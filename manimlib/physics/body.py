import numpy as np

from manimlib.constants import DIMENSIONS
from manimlib.mobject.mobject import Mobject
from manimlib.mobject.geometry import Polyline


class Body:
    """
    Represents a body in a physical system
    """

    def __init__(
        self,
        mass: float=1.0,
        charge: float=0.0,
        position: np.ndarray=np.array([0,0,0]),
        velocity: np.ndarray=np.array([0,0,0]),
        mobj: Mobject=None,
        tracer: Polyline=None
    ) -> None:
        """
        Initialize a new Body object

        Keyword arguments
        -----------------
        mass (float): the mass of the body (default: 1)
        position (np.ndarray[3]): vector representing the (x,y,z) components
                 of the position (default: [0,0,0])
        velocity (np.ndarray[3]): vector representing the (x,y,z) components
                 of the velocity (default: [0,0,0])
        mobj (Mobject): the Mobject this body is linked to (default: None)
        tracer (Polyline): polyline to trace the body's movement, must already
               have a vertex at the position of the body (default: None)
        """

        self.mass: float = mass
        if self.mass < 0:
            raise Exception(
                f"Current value for mass ({self.mass}) is negative!"
            )
        
        self.charge: float = charge
        
        self.position: np.ndarray = position
        if len(self.position.shape) != 1 or self.position.shape[0] != DIMENSIONS:
            raise Exception(f"Position {self.position} is invalid!")

        self.velocity: np.ndarray = velocity
        if len(self.velocity.shape) != 1 or self.velocity.shape[0] != DIMENSIONS:
            raise Exception(f"Velocity {self.velocity} is invalid!")
        
        self.mobj: Mobject = mobj
        if self.mobj is not None:
            self.mobj.move_to(self.position)
        
        if tracer is not None and not isinstance(tracer, Polyline):
            raise Exception(f"The tracer instance {tracer} provided is not a Polyline")
        self.tracer: Polyline = tracer

        # index of the body in the physical system, to be set
        # by a PhysicalSystem
        self.index: int = -1

        # Save the movement of the body during animation to update the tracer, due to a problem
        # with obtaining vertices from an existing polyline
        self.path: list[np.ndarray] = [self.position]

    def __str__(self) -> str:
        return (f"{self.__class__.__name__}<mass={self.mass},"
                f"position={self.position},velocity={self.velocity},"
                f"mobject={self.mobj},tracer={self.tracer},"
                f"index={self.index}>") 

    def set_mass(self, mass: float) -> None:
        """
        Set the mass of the body

        Keyword arguments
        -----------------
        mass (float): the new mass
        """
        self.mass = mass
        if self.mass < 0:
            raise Exception(
                f"Current value for mass ({self.mass}) is negative!"
            )
    
    def set_charge(self, charge: float) -> None:
        """
        Set the charge of the body

        Keyword arguments
        -----------------
        charge (float): the new charge
        """
        self.charge = charge

    def set_position(
        self,
        position: np.ndarray,
        update_mobject_position: bool=False
    ) -> None:
        """
        Set the position of the body (and move its mobject there)

        Keyword arguments
        -----------------
        position (np.ndarray[3]): vector representing the (x,y,z) components
                 of the position
        update_mobject_position (bool): whether to update the position of the
                                mobject to the new position (default: False)
        """
        self.position = position
        if len(self.position.shape) != 1 or self.position.shape[0] != DIMENSIONS:
            raise Exception(f"Position {self.position} is invalid!")
        if update_mobject_position:
            self.update_mobject_position()
    
    def update_mobject_position(self) -> None:
        """
        Update the mobject's position to the current body position
        """
        if self.mobj is not None:
            self.mobj.move_to(self.position)
    
    def update_tracer(self) -> None:
        """
        Update the tracer mobject for the body
        """
        if self.tracer is not None:
            self.path.append(self.position)
            self.tracer.become(
                Polyline(
                    *self.path,
                    stroke_color=self.tracer.get_stroke_color(),
                    stroke_opacity=self.tracer.get_stroke_opacity(),
                    stroke_width=self.tracer.get_stroke_width()
                )
            )
    
    def set_velocity(self, velocity: np.ndarray) -> None:
        """
        Set the velocity of the body

        Keyword arguments
        -----------------
        velocity (np.ndarray[3]): vector representing the (x,y,z) components
                 of the velocity
        """
        self.velocity = velocity
        if len(self.velocity.shape) != 1 or self.velocity.shape[0] != DIMENSIONS:
            raise Exception(f"Velocity {self.velocity} is invalid!")
    
    def set_mobj(self, mobj: Mobject) -> None:
        """
        Set a mobject for this body
        (and moves it to the body position)

        Keyword arguments
        -----------------
        mobj (Mobject): the Mobject instance
        """
        self.mobj = mobj
        if self.mobj is not None:
            self.mobj.move_to(self.position)

