import numpy as np

from manimlib.constants import TAU

def angle_to_full_circle(angle: float) -> float:
    """
    Transform angle [-pi, pi] to [0, 2pi] range

    Keyword arguments
    -----------------
    angle (float): the angle in radians to transform

    Returns
    -----------------
    The transformed angle (float)
    """
    return angle if angle >= 0 else TAU + angle

def arctan2(y: float, x: float) -> float:
    """
    Arctangent of angle (respects quadrants)

    Keyword arguments
    -----------------
    y (float): y component
    x (float): x component

    Returns
    -----------------
    the angle (float) in radians
    """
    return angle_to_full_circle(np.arctan2(y, x))
