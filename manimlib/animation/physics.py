from __future__ import annotations

import numpy as np

from typing import TYPE_CHECKING

from manimlib.mobject.mobject import Mobject
from manimlib.scene.scene import Scene
from manimlib.animation.animation import Animation
from manimlib.physics.physical_system import PhysicalSystem
from manimlib.physics.integrator import symplectic_euler
from manimlib.constants import DIMENSIONS, DEFAULT_FPS

if TYPE_CHECKING:
    from typing import Callable


class EvolvePhysicalSystem(Animation):
    """
    Evolves a PhysicalSystem object
    (integrates it over time)
    """

    def __init__(
        self,
        mobject: PhysicalSystem,
        integrator: Callable[
            [
                PhysicalSystem,
                np.ndarray
            ],
            tuple[np.ndarray, np.ndarray]
        ]=symplectic_euler,
        t: np.ndarray=np.linspace(0, 10, DEFAULT_FPS*100),
        verbose: int=1,
        scene: Scene=None,
        background_mobjects: tuple[Mobject]=(),
        foreground_mobjects: tuple[Mobject]=(),
        **kwargs
    ):
        """
        Initialize a new EvolvePhysicalSystem instance

        Keyword arguments
        -----------------
        mobject (PhysicalSystem): the PhysicalSystem mobject
        integrator (Callable[
                [
                    PhysicalSystem,
                    np.ndarray,
                    int
                ],
                tuple[np.ndarray, np.ndarray]
            ]): function in terms of (system, t, verbose) returning two tensors (pos, vel) where
                pos[i,j,k] represents the component along the k-th axis of the position of
                body j at time-point i (same holds por vel).
                For an example of the 'integrator' function, see symplectic_euler in
                manimlib.physics.integrator, which is the default value.
        t (np.ndarray[N]): vector containing all time-points of integration, must be monotonic
                           (default np.linspace(0, 10, DEFAULT_FPS*100))
        verbose (int): level of verbosity of the integrator (default: 1)
        scene (Scene): scene is needed if we want to keep the body mobject, force mobject,
                       body tracer rendering order (default: None)
        background_mobjects (tuple[Mobject]): mobjects that form part of the background and should be
                            brought to the back of the rendering order in each animation step, first in the
                            tuple is the furthest object in the back (default: empty tuple)
        foreground_mobjects (tuple[Mobject]): mobjects that form part of the foreground and should be
                            brought to the front of the rendering order in each animation step, last in the
                            tuple is the closest object in the front (default: empty tuple)
        kwargs (dict[str, Any]): arguments to be interpreted by
               the Animation superclass
        """
        if not isinstance(mobject, PhysicalSystem):
            raise Exception(
                f"({self.__class__.__name__}) The mobject object "
                "passed to the constructor is not an instance of "
                "PhysicalSystem"
            )
        super().__init__(mobject, **kwargs)
        # Save properties
        self.n_bodies: int = self.mobject.get_n_bodies()
        self.n_tpoints: int = t.shape[0]
        self.scene: Scene = scene
        self.background_mobjects: tuple[Mobject] = background_mobjects
        self.foreground_mobjects: tuple[Mobject] = foreground_mobjects
        # Integrate the system
        self.positions, self.velocities = integrator(self.mobject, t, verbose)

    def interpolate_mobject(self, alpha: float) -> None:
        row_index: int = min(
            int(np.floor(alpha*self.n_tpoints)),
            self.n_tpoints-1
        )
        self.mobject.set(
            positions=self.positions[row_index],
            velocities=self.velocities[row_index]
        )
        # Update mobjects in the system
        self.mobject.update_mobjects(
            self.scene,
            self.background_mobjects,
            self.foreground_mobjects
        )
