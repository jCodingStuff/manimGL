from __future__ import annotations

import numpy as np
from manimlib.constants import DIMENSIONS

from manimlib.physics.physical_system import PhysicalSystem


def symplectic_euler(system: PhysicalSystem, t: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    """
    Perform symplectic euler integration of the physical system

    Keyword arguments
    -----------------
    system (PhysicalSystem): the system to integrate
    t (np.ndarray[N]): vector containing all time-points of integration, must be monotonic

    Returns
    -----------------
    positions (np.ndarray[N, n_bodies, DIMENSIONS]): the position of each body for every
              time-point in t
    velocities (np.ndarray[N, n_bodies, DIMENSIONS]): the velocity of each body for every
              time-point in t
    """
    n_bodies: int = system.get_n_bodies()
    n_tpoints: int = t.shape[0]
    positions: np.ndarray = np.zeros((n_tpoints, n_bodies, DIMENSIONS))
    velocities: np.ndarray = np.zeros((n_tpoints, n_bodies, DIMENSIONS))
    # Set initial positions and velocities
    positions[0] = system.get_positions()
    velocities[0] = system.get_velocities()
    # Integrate
    for i in range(1, n_tpoints):
        # Compute dt (delta time)
        dt: float = t[i] - t[i-1]
        # Compute accelerations
        system.update_positions(positions[i-1])
        system.update_velocities(velocities[i-1])
        accelerations: np.ndarray = system.compute_accelerations()
        # Record new velocities
        velocities[i] = velocities[i-1] + dt * accelerations
        # Record new positions
        positions[i] = positions[i-1] + dt * velocities[i]
    # Set the initial state back in the system
    system.update_positions(positions[0])
    system.update_velocities(velocities[0])
    # Return
    return positions, velocities