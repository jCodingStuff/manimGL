from __future__ import annotations

from manimlib.constants import DIMENSIONS
from manimlib.scene.scene import Scene
from manimlib.mobject.mobject import Group
from manimlib.physics.body import Body
from manimlib.physics.force import Force, NewtonGravitationalForce
import numpy as np

from typing import TYPE_CHECKING

if TYPE_CHECKING:
    from manimlib.mobject.mobject import Mobject


class PhysicalSystem(Group):
    """
    Represents a physical system
    (all submobjects in this class should have already been added
    to the scene)
    """

    def __init__(self, bodies: list[Body]=[], forces: list[Force]=[], **kwargs):
        """
        Initialize a new PhysicalSystem object

        Keyword arguments
        -----------------
        bodies (list[Body]): bodies to be added to the system (default: empty list)
        forces (list[Force]): forces to be added to the system (default: empty list)
        kwargs (dict[str, Any]): arguments to pass to the Group superclass
        """
        # Aggregate all submobjects and pass them to the superclass constructor
        mobjects: list[Mobject] = []
        for body in bodies:
            if body.mobj is not None:
                mobjects.append(body.mobj)
            if body.tracer is not None:
                mobjects.append(body.tracer)
        for force in forces:
            for mobj in force.mobjects:
                if mobj is not None:
                    mobjects.append(mobj)
        super().__init__(*mobjects, **kwargs)
        self.bodies: list[Body] = bodies
        for i, body in enumerate(self.bodies):  # set indices for the bodies
            body.index = i
        self.forces: list[Force] = forces

    def __str__(self):
        body_str: str = [str(body) for body in self.bodies]
        force_str: str = [str(force) for force in self.forces]
        return f"{self.__class__.__name__}<bodies={body_str},forces={force_str}>"

    def update_mobjects(
        self,
        scene: Scene,
        background_mobjects: tuple[Mobject]=(),
        foreground_mobjects: tuple[Mobject]=()
    ) -> None:
        """
        Update the mobjects in the system so that they
        are correctly displayed

        Keyword arguments
        -----------------
        scene (Scene): the scene in which the mobjects are
        background_mobjects (tuple[Mobject]): mobjects that form part of the background and should be
                            brought to the back of the rendering order in each animation step, first in the
                            tuple is the furthest object in the back (default: empty tuple)
        foreground_mobjects (tuple[Mobject]): mobjects that form part of the foreground and should be
                            brought to the front of the rendering order in each animation step, last in the
                            tuple is the closest object in the front (default: empty tuple)
        """
        for body in self.bodies:
            body.update_tracer()
        # Update mobjects in forces
        for force in self.forces:
            force.update_mobjects()
        # Update body mobject
        for body in self.bodies:
            body.update_mobject_position()
            if body.mobj is not None:
                scene.bring_to_front(body.mobj)
        # Background mobjects
        if background_mobjects:
            scene.bring_to_back(*background_mobjects)
        # Foreground mobjects
        if foreground_mobjects:
            scene.bring_to_front(*foreground_mobjects)

    def get_n_bodies(self) -> int:
        """
        Get the number of bodies in the system

        Returns
        -----------------
        The number of bodies (int)
        """
        return len(self.bodies)

    def get_n_forces(self) -> int:
        """
        Get the number of forces in the system

        Returns
        -----------------
        The number of forces (int)
        """
        return len(self.forces)

    def get_masses(self) -> np.ndarray:
        """
        Get a vector containing the masses of the bodies in
        the system

        Returns
        -----------------
        An np.ndarray[n_bodies] of floats containing the masses
        """
        return np.array([body.mass for body in self.bodies])

    def compute_forces(self) -> np.ndarray:
        """
        Compute the overall force exerted on each body in the system

        Returns
        -----------------
        An np.ndarray[n_bodies, DIMENSIONS] of floats containing
        the forces
        """
        forces: np.ndarray = np.zeros((self.get_n_bodies(), DIMENSIONS))
        for force in self.forces:
            force.apply(forces)
        return forces
    
    def compute_accelerations(self) -> np.ndarray:
        """
        Compute the overall acceleration resulting on each body in
        the system

        Returns
        -----------------
        An np.ndarray[n_bodies, DIMENSIONS] of floats containing
        the accelerations
        """
        return self.compute_forces() / self.get_masses().reshape((self.get_n_bodies(), 1))

    def get_positions(self) -> np.ndarray:
        """
        Get the position of each body

        Returns
        -----------------
        An np.ndarray[n_bodies, DIMENSIONS] with the positions
        """
        positions: np.ndarray = np.zeros((self.get_n_bodies(), DIMENSIONS))
        for i, body in enumerate(self.bodies):
            positions[i] = body.position
        return positions

    def set_positions(
        self,
        positions: np.ndarray,
        update_mobjects_position: bool=False
    ) -> None:
        """
        Update positions of the bodies in the system

        Keyword arguments
        -----------------
        positions (np.ndarray[n_bodies, DIMENSIONS]): position for
                  each body in the system
        update_mobjects_position (bool): True if we want to update the position
                                 of the mobject of each body, False otherwise
                                 (default: False)
        """
        desired_shape: tuple[int, int] = (self.get_n_bodies(), DIMENSIONS)
        if positions.shape != desired_shape:
            raise Exception(
                f"The shape of the provided positions {positions.shape}"
                f" does not match {desired_shape}"
            )
        for body, position in zip(self.bodies, positions):
            body.set_position(position, update_mobject_position=update_mobjects_position)
    
    def get_velocities(self) -> np.ndarray:
        """
        Get the velocity of each body

        Returns
        -----------------
        An np.ndarray[n_bodies, DIMENSIONS] with the velocities
        """
        velocities: np.ndarray = np.zeros((self.get_n_bodies(), DIMENSIONS))
        for i, body in enumerate(self.bodies):
            velocities[i] = body.velocity
        return velocities

    def set_velocities(self, velocities: np.ndarray) -> None:
        """
        Update velocities of the bodies in the system

        Keyword arguments
        -----------------
        velocities (np.ndarray[n_bodies, DIMENSIONS]): position for
                  each body in the system
        """
        desired_shape: tuple[int, int] = (self.get_n_bodies(), DIMENSIONS)
        if velocities.shape != desired_shape:
            raise Exception(
                f"The shape of the provided velocities {velocities.shape}"
                f" does not match {desired_shape}"
            )
        for body, velocity in zip(self.bodies, velocities):
            body.set_velocity(velocity)

    def set(
        self,
        **kwargs
    ) -> None:
        """
        Set properties of the system

        Keyword arguments
        -----------------
        kwargs (dict[str, Any]): the key determines the name of the
               property to set, whose existance is checked by ensuring
               that there exists an instance method with the name
               'set_<key>'
        """
        for k, v in kwargs.items():
            setter_name: str = f"set_{k}"
            if hasattr(self, setter_name):
                setter = getattr(self, setter_name)
                if callable(setter):
                    setter(v)
                else:
                    raise Exception(
                        f"({self.__class__.__name__}) {setter_name} "
                        f"found but not callable"
                    )
            else:
                raise Exception(
                    f"({self.__class__.__name__}) No instance "
                    f"member named {setter_name} found"
                )

    def fill_forces(self, **kwargs) -> None:
        """
        Some subclasses will be able to fill forces by themselves with 
        this method

        Keyword arguments
        -----------------
        kwargs (dict[str, Any]): option dictionary to be interpreted by each
               subclass
        """
        pass


class GravitationalSystem(PhysicalSystem):
    """
    Special case of a physical system for gravitational simulation
    """

    def fill_forces(self, **kwargs) -> None:
        """
        Add a NewtonGravitationalForce to every pair of bodies
        TODO: add an option in 'kwargs' to automatically add line mobjects for the force
              (they will have to be added with the self.add() method as well, since this is
              a mobject itself too)

        Keyword arguments
        -----------------
        kwargs (dict[str, Any]): understood keys are 'model' (default: 'newton',
                and others in the future) and 'G' (default: 1.0, gravitational
                constant)
        """
        # Check if there are at least 2 bodies
        if self.get_n_bodies() <= 1:
            return
        # Set default value for 'model'
        if "model" not in kwargs:
            kwargs["model"] = 'newton'
        # Add the forces
        if kwargs["model"] == "newton":
            # Get the G value
            if 'G' not in kwargs:
                kwargs['G'] = 1.0
            for i, body1 in enumerate(self.bodies):
                for _, body2 in enumerate(self.bodies[i+1:],start=i+1):
                    self.forces.append(
                        NewtonGravitationalForce(
                            (body1, body2),
                            G = kwargs["G"]
                        )
                    )
        else:
            raise Exception(f"Force model with name {kwargs['model']} could not be found...")