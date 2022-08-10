from __future__ import annotations

import numpy as np
from manimlib.constants import DIMENSIONS

from manimlib.physics.physical_system import PhysicalSystem


def symplectic_euler(
    system: PhysicalSystem,
    t: np.ndarray,
    verbose: int=0
) -> tuple[np.ndarray, np.ndarray]:
    """
    Perform symplectic euler integration of the physical system

    Keyword arguments
    -----------------
    system (PhysicalSystem): the system to integrate
    t (np.ndarray[N]): vector containing all time-points of integration, must be monotonic
    verbose (int): level of verbosity (default: 0)

    Returns
    -----------------
    positions (np.ndarray[N, n_bodies, DIMENSIONS]): the position of each body for every
              time-point in t
    velocities (np.ndarray[N, n_bodies, DIMENSIONS]): the velocity of each body for every
              time-point in t
    """
    # For verbose stuff
    percentages: list[int] = []
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
        system.set(
            positions=positions[i-1],
            velocities=velocities[i-1]
        )
        accelerations: np.ndarray = system.compute_accelerations()
        # Record new velocities
        velocities[i] = velocities[i-1] + dt * accelerations
        # Record new positions
        positions[i] = positions[i-1] + dt * velocities[i]
        # Print stuff if verbose
        if verbose:
            percentage: int = round(i/(n_tpoints-1)*100)
            if percentage not in percentages and percentage%25==0:
                percentages.append(percentage)
                print(f"Integration progress: {percentage}%")
    # Set the initial state back in the system
    system.set(
        positions=positions[0],
        velocities=velocities[0]
    )
    # Return
    return positions, velocities

def rk4(system: PhysicalSystem, t: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    """
    Perform 4th order Runge-Kutta integration of the physical system

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
        # k1
        # k2
        # k3
        # k4
        # Aggregate all components
    # Set the initial state back in the system
    system.set_positions(positions[0])
    system.set_velocities(velocities[0])
    # Return
    return None, None