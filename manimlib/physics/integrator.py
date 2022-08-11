from __future__ import annotations

import numpy as np
from manimlib.constants import DIMENSIONS

from manimlib.physics.physical_system import PhysicalSystem
from manimlib.utils.simple_functions import printProgressBar


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
    if verbose:
        print("Integrating with Symplectic Euler...")

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
            printProgressBar(i, n_tpoints-1, "Progress: ", "Complete", length=50)
    
    # Set the initial state back in the system
    system.set(
        positions=positions[0],
        velocities=velocities[0]
    )
    # Return
    return positions, velocities

def rk4(
    system: PhysicalSystem,
    t: np.ndarray,
    verbose: int=0
) -> tuple[np.ndarray, np.ndarray]:
    """
    Perform 4th order Runge-Kutta integration of the physical system
    Not adequate for stochastic systems. This is due to the fact that RK4 computes
    the new state with a weighted average of "present and future" systems, thus
    smoothing the effect of stochastic forces.

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
    if verbose:
        print("Integrating with 4th order Runge-Kutta...")

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
        system.set(
            positions=positions[i-1],
            velocities=velocities[i-1]
        )
        k1 = (velocities[i-1], system.compute_accelerations())
        # k2
        k2_vel = velocities[i-1]+dt*k1[1]/2
        system.set(
            positions=positions[i-1]+dt*k1[0]/2,
            velocities=k2_vel
        )
        k2 = (k2_vel, system.compute_accelerations())
        # k3
        k3_vel = velocities[i-1]+dt*k2[1]/2
        system.set(
            positions=positions[i-1]+dt*k2[0]/2,
            velocities=k3_vel
        )
        k3 = (k3_vel, system.compute_accelerations())
        # k4
        k4_vel = velocities[i-1]+dt*k3[1]
        system.set(
            positions=positions[i-1]+dt*k3[0],
            velocities=k4_vel
        )
        k4 = (k4_vel, system.compute_accelerations())
        # Aggregate all components and store new positions and velocities
        positions[i] = positions[i-1] + dt*(k1[0]+2*k2[0]+2*k3[0]+k4[0])/6
        velocities[i] = velocities[i-1] + dt*(k1[1]+2*k2[1]+2*k3[1]+k4[1])/6
        # Print stuff if verbose
        if verbose:
            printProgressBar(i, n_tpoints-1, "Progress: ", "Complete", length=50)
    
    # Set the initial state back in the system
    system.set(
        positions=positions[0],
        velocities=velocities[0]
    )
    # Return
    return positions, velocities
