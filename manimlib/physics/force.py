from __future__ import annotations

from abc import ABCMeta, abstractmethod

import numpy as np

from typing import Union
from manimlib.constants import DIMENSIONS, KB, PI, ROOM_TEMPERATURE, TAU, EPS0, NO_SEED
from manimlib.mobject.mobject import Mobject
from manimlib.utils.math import arctan2

from manimlib.physics.body import Body
from manimlib.mobject.geometry import Line, Arc
from manimlib.mobject.three_dimensions import Line3D


class Force(metaclass=ABCMeta):
    """
    Represents a force among/for bodies in a physical
    system
    """

    def __init__(
        self,
        bodies: tuple[Body, ...],
        mobjects: tuple[Mobject, ...]=()
    ) -> None:
        """
        Initialize a new Force object

        Keyword arguments
        -----------------
        bodies (tuple[Body, ...]): the bodies to which the force applies
        mobjects (tuple[Mobject, ...]): mobject(s) representing the force, should
                 be already set to a desired position and added to the scene
                 (subclasses know how to update them). Regular shapes should
                 be used for 2D simulations and 3D shapes for 3D simulations!
                 Otherwise things MAY BREAK! (default: empty tuple)
        """
        self.bodies: tuple[Body, ...] = bodies
        if not self.bodies:
            raise Exception("No bodies have been provided!")
        self.mobjects: tuple[Mobject, ...] = mobjects
    
    def __str__(self) -> str:
        body_info: str = ','.join([f"body{i}_index={body.index}" for i, body in enumerate(self.bodies, start=1)])
        return (f"{body_info},mobjects={self.mobjects}")

    def update_mobjects(self) -> None:
        """
        Update the foce mobject(s)
        (to be implemented by each subclass)
        """
        pass

    @abstractmethod
    def apply(self, forces: np.ndarray) -> None:
        """
        Apply the force and add the contributions to the 'total'
        forces in the system

        Keyword arguments
        -----------------
        forces (np.ndarray[nbodies, DIMENSIONS]): each row stores the
               total force exerted in the body with that index
        """
        pass


class SingleForce(Force):
    """
    Abstract force class (apply() is not implemented) that applies to
    1 body
    """
    def __init__(
        self,
        bodies: tuple[Body],
        mobjects: tuple[Mobject, ...]=()
    ) -> None:
        """
        Initialize a new SingleForce object

        Keyword arguments
        -----------------
        bodies (tuple[Body]): the body to which the force applies
        mobjects (tuple[Mobject, ...]): mobject(s) representing the force, should
                 be already set to a desired position (subclasses know
                 how to update them). Regular shapes should be used for 2D
                 simulations and 3D shapes for 3D simulations! Otherwise
                 things MAY BREAK! (default: empty tuple)
        """
        if len(bodies) != 1:
            raise Exception("You must provide exactly 1 body!")
        super().__init__(bodies, mobjects)


class PairForce(Force):
    """
    Abstract force class (apply() is not implemented) that applies to
    2 bodies
    """
    def __init__(
        self,
        bodies: tuple[Body, Body],
        mobjects: tuple[Mobject, ...]=(),
    ) -> None:
        """
        Initialize a new PairForce object

        Keyword arguments
        -----------------
        bodies (tuple[Body, Body]): the 2 bodies to which the force applies
        mobjects (tuple[Mobject, ...]): mobject(s) representing the force, should
                 be already set to a desired position (subclasses know
                 how to update them). Regular shapes should be used for 2D
                 simulations and 3D shapes for 3D simulations! Otherwise
                 things MAY BREAK! (default: empty tuple)
        """
        if len(bodies) != 2:
            raise Exception("You must provide exactly 2 bodies!")
        super().__init__(bodies, mobjects)


class TripletForce(Force):
    """
    Abstract force class (apply() is not implemented) that applies to
    2 bodies
    """
    def __init__(
        self,
        bodies: tuple[Body, Body, Body],
        mobjects: tuple[Mobject, ...]=(),
    ) -> None:
        """
        Initialize a new PairForce object

        Keyword arguments
        -----------------
        bodies (tuple[Body, Body, Body]): the 3 bodies to which the force applies
        mobjects (tuple[Mobject, ...]): mobject(s) representing the force, should
                 be already set to a desired position (subclasses know
                 how to update them). Regular shapes should be used for 2D
                 simulations and 3D shapes for 3D simulations! Otherwise
                 things MAY BREAK! (default: empty tuple)
        """
        if len(bodies) != 3:
            raise Exception("You must provide exactly 3 bodies!")
        super().__init__(bodies, mobjects)


class PairLineForce(PairForce):
    """
    Abstract force class (apply() is not implemented) that can manage a
    line connecting its 2 bodies as mobject
    """

    def __init__(
        self,
        bodies: tuple[Body, Body],
        mobjects: tuple[Union[Line, Line3D]]=(),
    ) -> None:
        """
        Initialize a new PairLineForce instance

        Keyword arguments
        -----------------
        bodies (tuple[Body, Body]): the 2 bodies to which the force applies
        mobjects (tuple[Line | Line3D]): line representing the force and follows
                 the bodies, may already be set to have the bodies at its
                 extremes (start=bodies[0], end=bodies[1]) (default: empty tuple)
        """
        if mobjects and not isinstance(mobjects[0], (Line, Line3D)):
            raise Exception(f"The provided mobject {mobjects[0]} is not a Line/Line3D")
        super().__init__(bodies, mobjects)
    
    def update_mobjects(self) -> None:
        """
        Update the line mobject to follow the bodies
        """
        if not self.mobjects:
            return
        line: Union[Line, Line3D] = self.mobjects[0]
        body1, body2 = self.bodies
        if isinstance(line, Line):  # for 2D simulations
            line.become(
                Line(
                    body1.position,
                    body2.position,
                    stroke_color=line.get_stroke_color(),
                    stroke_opacity=line.get_stroke_opacity(),
                    stroke_width=line.get_stroke_width()
                )
            )
        elif isinstance(line, Line3D):  # TODO: TEST THIS
            line.become(
                Line3D(
                    body1.position,
                    body2.position,
                    width=line.get_width(),
                    color=line.get_color(),
                    opacity=line.get_opacity()
                )
            )
        else:  # This should never happen
            print(
                f"({self.__class__.__name__}) WARNING: the mobject "
                "is not a line. Cannot update it."
            )


class TripletArcForce(TripletForce):
    """
    Abstract force class (apply() is not implemented) that can manage an
    arc representing the body1-body2-body3 angle
    """
    def __init__(
        self,
        bodies: tuple[Body, Body, Body],
        mobjects: tuple[Arc]=(),
    ) -> None:
        """
        Initialize a new TripletArcForce instance

        Keyword arguments
        -----------------
        bodies (tuple[Body, Body, Body]): the 3 bodies to which the force applies
        mobjects (tuple[Arc]): arc representing the force and defines the angle
                 body1-body2-body3, may already be set to have its desired starting
                 position) (default: empty tuple)
        """
        if mobjects and not isinstance(mobjects[0], Arc):
            raise Exception(f"The provided mobject {mobjects[0]} is not an Arc")
        super().__init__(bodies, mobjects)
    
    def update_mobjects(self) -> None:
        """
        Update the arc to follow the angle
        """
        if not self.mobjects:
            return
        arc: Arc = self.mobjects[0]
        body1, body2, body3 = self.bodies
        r12: np.ndarray = body1.position - body2.position
        r32: np.ndarray = body3.position - body2.position
        angle12: float = arctan2(r12[1], r12[0])
        angle32: float = arctan2(r32[1], r32[0])
        start_angle: float = min(angle12, angle32)
        end_angle: float = max(angle12, angle32)
        angle: float = end_angle - start_angle
        if angle > PI:  # Reverse > PI angles
            angle = angle - TAU
        arc.become(
            Arc(
                start_angle,
                angle,
                radius=arc.radius,
                arc_center=body2.position,
                stroke_color=arc.get_stroke_color(),
                stroke_opacity=arc.get_stroke_opacity(),
                stroke_width=arc.get_stroke_width(),
                fill_color=arc.get_fill_color(),
                fill_opacity=arc.get_fill_opacity()
            )
        )


class LangevinFrictionForce(SingleForce):
    """
    Friction force according to Langevin dynamics:
    m_i * a_i = F_i = ... - gamma * m_i * v_i,
    where gamma is the friction coefficient
    """
    def __init__(
        self,
        bodies: tuple[Body],
        mobjects: tuple[Mobject, ...]=(),
        gamma: float=1.0
    ) -> None:
        """
        Initialize a new LangevinFrictionForce object

        Keyword arguments
        -----------------
        bodies (tuple[Body]): the body to which the force applies
        mobjects (tuple[Mobject, ...]): mobject(s) representing the force, should
                 be already set to a desired position (subclasses know
                 how to update them). Regular shapes should be used for 2D
                 simulations and 3D shapes for 3D simulations! Otherwise
                 things MAY BREAK! (default: empty tuple)
        gamma (float): friction coefficient (default: 1.0)
        """
        self.gamma: float = gamma
        super().__init__(bodies, mobjects)
    
    def __str__(self) -> str:
        return (f"{self.__class__.__name__}<{super().__str__()},"
                f"gamma={self.gamma}>")
    
    def apply(self, forces: np.ndarray) -> None:
        body: Body = self.bodies[0]
        forces[body.index] -= self.gamma * body.mass * body.velocity


class LangevinHeatBathForce(SingleForce):
    """
    Force simulating the body colliding with particles in a heat bath
    (at a certain temperature)
    According to the Langevin equation of motion:
    m_i * a_i = F_i = ... + R_i, where R_i is a vector where the components
    come from a normal distribution with mean zero and variance
    2 * m_i * gamma * kb * T. gamma is the friction coefficient, kb is
    the Boltzmann constant, and T is the temperature
    """

    def __init__(
        self,
        bodies: tuple[Body],
        mobjects: tuple[Mobject, ...]=(),
        tridimensional: bool=False,
        seed: int=NO_SEED,
        gamma: float=1.0,
        kb: float=KB,
        T: float=ROOM_TEMPERATURE
    ) -> None:
        """
        Initialize a new LangevinHeatBathForce object

        Keyword arguments
        -----------------
        bodies (tuple[Body]): the body to which the force applies
        mobjects (tuple[Mobject, ...]): mobject(s) representing the force, should
                 be already set to a desired position (subclasses know
                 how to update them). Regular shapes should be used for 2D
                 simulations and 3D shapes for 3D simulations! Otherwise
                 things MAY BREAK! (default: empty tuple)
        tridimensional (bool): False if system is 2D, True if 3D (default: False)
        seed (int): seed for the random number generator. If -1, then no seed is
             passed to the generator, yielding non-reproducible results (default: -1)
        gamma (float): friction coefficient (default: 1.0)
        kb (float): Boltzmann constant (default: 1.380649e-23)
        T (float): temperature (default: 293.0)
        """
        self.tridimensional: bool = tridimensional
        self.generator: np.random.Generator = np.random.default_rng(seed=None if seed==NO_SEED else seed)
        self.gamma: float = gamma
        self.kb: float = kb
        self.T: float = T
        super().__init__(bodies, mobjects)
    
    def __str__(self) -> str:
        return (f"{self.__class__.__name__}<{super().__str__()},"
                f"tridimensional={self.tridimensional},"
                f"generator={self.generator},gamma={self.gamma},"
                f"kb={self.kb},T={self.T}>")
    
    def apply(self, forces: np.ndarray) -> None:
        body: Body = self.bodies[0]
        # We work out the deviation every time just in case in the future
        # we want to perform simulations with variable temperature
        deviation: float = np.sqrt(2*body.mass*self.gamma*self.kb*self.T)
        force: np.ndarray = self.generator.normal(scale=deviation, size=DIMENSIONS)
        if not self.tridimensional:
            force[2] == 0.0
        forces[body.index] += force
        

class NewtonGravitationalForce(PairLineForce):
    """
    Newton's gravitational force between two masses
    F = - G * m1 * m2 / r^2
    """

    def __init__(
        self,
        bodies: tuple[Body, Body],
        mobjects: tuple[Union[Line, Line3D]]=(),
        G: float=1.0
    ) -> None:
        """
        Initialize a new NewtonGravitationalForce instance

        Keyword arguments
        -----------------
        bodies (tuple[Body, Body]): the bodies to which the force applies
        mobjects (tuple[Line | Line3D]): line representing the force and follows
             the bodies, must already be set to have the bodies at its
             extremes (start=bodies[0], end=bodies[1]) (default: empty tuple)
        G (float): gravitational constant (default: 1.0)
        """
        super().__init__(bodies, mobjects)
        self.G: float = G
    
    def __str__(self) -> str:
        return (f"{self.__class__.__name__}<{super().__str__()},"
                f"G={self.G}>")

    def apply(self, forces: np.ndarray) -> None:
        body1, body2 = self.bodies
        delta_pos: np.ndarray = body2.position - body1.position
        distance: float = np.linalg.norm(delta_pos)
        if distance > 0:  # apply the force
            force: np.ndarray = self.G * body1.mass * body2.mass * delta_pos / distance**3
            forces[body1.index] += force
            forces[body2.index] -= force


class HarmonicBondForce(PairLineForce):
    """
    Harmonic bond force between two bonded atoms
    """

    def __init__(
        self,
        bodies: tuple[Body, Body],
        mobjects: tuple[Union[Line, Line3D]]=(),
        k: float=1.0,
        r0: float=1.0
    ) -> None:
        """
        Initialize a new HarmonicBondForce instance

        Keyword arguments
        -----------------
        bodies (tuple[Body, Body]): the bodies to which the force applies
        mobjects (tuple[Line | Line3D]): line representing the force and follows
             the bodies, must already be set to have the bodies at its
             extremes (start=bodies[0], end=bodies[1]) (default: empty tuple)
        k (float): force constant (default: 1.0)
        r0 (float): equilibrium length (default: 1.0)
        """
        super().__init__(bodies, mobjects)
        self.k: float = k
        self.r0: float = r0
    
    def __str__(self) -> str:
        return (f"{self.__class__.__name__}<{super().__str__()},"
                f"k={self.k},r0={self.r0}>")
    
    def apply(self, forces: np.ndarray) -> None:
        body1, body2 = self.bodies
        r12: np.ndarray = body1.position - body2.position
        dist12: float = np.linalg.norm(r12)
        if dist12 > 0:
            force: np.ndarray = self.k * (dist12 - self.r0) * r12 / dist12
            forces[body1.index] -= force
            forces[body2.index] += force


class MorseBondForce(PairLineForce):
    """
    Morse bond force between two bonded atoms
    """

    def __init__(
        self,
        bodies: tuple[Body, Body],
        mobjects: tuple[Union[Line, Line3D]]=(),
        D: float=1.0,
        beta: float=1.0,
        r0: float=1.0
    ) -> None:
        """
        Initialize a new MorseBondForce instance

        Keyword arguments
        -----------------
        bodies (tuple[Body, Body]): the bodies to which the force applies
        mobjects (tuple[Line | Line3D]): line representing the force and follows
             the bodies, must already be set to have the bodies at its
             extremes (start=bodies[0], end=bodies[1]) (default: empty tuple)
        D (float): bond dissociation energy (default: 1.0)
        beta (float): force constant (default: 1.0)
        r0 (float): equilibrium length (default: 1.0)
        """
        super().__init__(bodies, mobjects)
        self.D: float = D
        self.beta: float = beta
        self.r0: float = r0
    
    def __str__(self) -> str:
        return (f"{self.__class__.__name__}<{super().__str__()},"
                f"D={self.D},beta={self.beta},r0={self.r0}>")

    def apply(self, forces: np.ndarray) -> None:
        body1, body2 = self.bodies
        r12: np.ndarray = body1.position - body2.position
        dist12: float = np.linalg.norm(r12)
        if dist12 > 0:
            exp_term: float = np.exp(-self.beta * (dist12 - self.r0))
            force: np.ndarray = 2 * self.D * self.beta * (1 - exp_term) * exp_term * r12 / dist12
            forces[body1.index] -= force
            forces[body2.index] += force


class HarmonicAngleForce(TripletArcForce):
    """
    Harmonic angle force among atoms 1,2,3
    with bonds 1-2 and 2-3
    """

    def __init__(
        self,
        bodies: tuple[Body, Body, Body],
        mobjects: tuple[Arc]=(),
        k: float=1.0,
        theta0: float=PI/2
    ) -> None:
        """
        Initialize a new HarmonicAngleForce instance

        Keyword arguments
        -----------------
        bodies (tuple[Body, Body, Body]): the bodies to which the force applies
        mobjects (tuple[Arc]): arc representing the force and defines the angle
                 body1-body2-body3, may already be set to have its desired starting
                 position) (default: empty tuple)
        k (float): force constant (default: 1.0)
        theta0 (float): equilibrium angle (default: PI/2)
        """
        super().__init__(bodies, mobjects)
        self.k: float = k
        self.theta0: float = theta0
    
    def __str__(self) -> str:
        return (f"{self.__class__.__name__}<{super().__str__()},"
                f"k={self.k},theta0={self.theta0}>")
    
    @staticmethod
    def are_collinear(
        point1: np.ndarray,
        point2: np.ndarray,
        point3: np.ndarray,
        tol: float=1e-6
    ) -> bool:
        """
        Check if 3 points are collinear

        Keyword arguments
        -----------------
        point1 (np.ndarray[3]): x,y,z coordinates of point1
        point2 (np.ndarray[3]): x,y,z coordinates of point2
        point3 (np.ndarray[3]): x,y,z coordinates of point3
        tol (float): tolerance (default: 1e-6)

        Returns
        -----------------
        True if they are collinear, False otherwise
        """
        dists: list[float] = [
            np.linalg.norm(point1-point2),
            np.linalg.norm(point1-point3),
            np.linalg.norm(point2-point3),
        ]
        dists.sort(reverse=True)
        return np.abs(dists[0] - np.sum(dists[1:])) <= tol
 
    def apply(self, forces: np.ndarray) -> None:
        body1, body2, body3 = self.bodies
        # Do not apply if angle is PI (180 degrees)
        if self.are_collinear(body1.position, body2.position, body3.position):
            return
        # Get angle
        r12: np.ndarray = body1.position - body2.position
        dist12: float = np.linalg.norm(r12)
        r32: np.ndarray = body3.position - body2.position
        dist32: float = np.linalg.norm(r32)
        unit_vector_1: np.ndarray = r12 / dist12
        unit_vector_2: np.ndarray = r32 / dist32
        cosine: float = np.dot(unit_vector_1, unit_vector_2)
        theta: float = np.arccos(cosine)
        # Common force terms
        common_term: float = self.k * (theta - self.theta0)
        inverse_term: float = 1/np.sqrt(1-cosine*cosine)
        # Body1
        force1: np.ndarray = -common_term * inverse_term * (r32/(dist12*dist32)-cosine*r12/(dist12*dist12))
        forces[body1.index] -= force1
        # Body3
        force3: np.ndarray = -common_term * inverse_term * (r12/(dist12*dist32)-cosine*r32/(dist32*dist32))
        forces[body3.index] -= force3
        # Body2
        forces[body2.index] += force1 + force3


class CoulombForce(PairForce):
    """
    Coulomb force between two bodies (1, 2)
    F1 = - f * (q1 * q2 / (eps_r * dist12^2)) * r12 / dist12^2
    f = 1 / (1*PI*eps_0)
    F2 = -F1

    eps_0 is the vacuum permittivity
    eps_r is the relative permittivity -> eps / esp_0, eps being the
                                          absolute permittivity of the
                                          dielectric
    q1 and q2 are the charges of bodies 1 and 2, respectively
    """
    def __init__(
        self,
        bodies: tuple[Body, Body],
        mobjects: tuple[Mobject, ...]=(),
        f: float=1/(4*PI*EPS0),
        eps_r: float=1.0
    ) -> None:
        """
        Initialize a new PairForce object

        Keyword arguments
        -----------------
        bodies (tuple[Body, Body]): the 2 bodies to which the force applies
        mobjects (tuple[Mobject, ...]): mobject(s) representing the force,
                 this subclass doesn't manage them (default: empty tuple)
        f (float): force constant together with eps_r (default: 1/(4*PI*EPS0))
        eps_r (float): relative permittivity, force constant together with f
              (default: 1.0) 
        """
        super().__init__(bodies, mobjects)
        self.f: float = f
        self.eps_r: float = eps_r
    
    def __str__(self) -> str:
        return (f"{self.__class__.__name__}<{super().__str__()},"
                f"f={self.f},eps_r={self.eps_r}>")
    
    def apply(self, forces: np.ndarray) -> None:
        body1, body2 = self.bodies
        r12: np.ndarray = body1.position - body2.position
        dist12: float = np.linalg.norm(r12)
        # Body1
        force1: np.ndarray = self.f * body1.charge * body2.charge * r12 / (self.eps_r * dist12**3)
        forces[body1.index] -= force1
        # Body2
        forces[body2.index] += force1
