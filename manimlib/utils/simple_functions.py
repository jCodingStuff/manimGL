from functools import lru_cache
import inspect
import math

import numpy as np


def sigmoid(x):
    return 1.0 / (1 + np.exp(-x))


@lru_cache(maxsize=10)
def choose(n, k):
    return math.comb(n, k)


def gen_choose(n, r):
    return np.prod(np.arange(n, n - r, -1)) / math.factorial(r)


def get_num_args(function):
    return len(get_parameters(function))


def get_parameters(function):
    return inspect.signature(function).parameters

# Just to have a less heavyweight name for this extremely common operation
#
# We may wish to have more fine-grained control over division by zero behavior
# in the future (separate specifiable values for 0/0 and x/0 with x != 0),
# but for now, we just allow the option to handle indeterminate 0/0.


def clip(a, min_a, max_a):
    if a < min_a:
        return min_a
    elif a > max_a:
        return max_a
    return a


def cap_magnitude(v: np.ndarray, max_magnitude: float) -> bool:
    """
    Cap (in-place) the magnitude of a vector by a given maximum value

    Keyword arguments
    -----------------
    v (np.ndarray): the vector
    max_magnitude (float): the allowed maximum magnitude

    Returns
    -----------------
    True if the vector was capped, false otherwise
    """
    magnitude: float = np.linalg.norm(v)
    if magnitude > max_magnitude:
        v *= max_magnitude/magnitude
        return True
    return False


def fdiv(a, b, zero_over_zero_value=None):
    if zero_over_zero_value is not None:
        out = np.full_like(a, zero_over_zero_value)
        where = np.logical_or(a != 0, b != 0)
    else:
        out = None
        where = True

    return np.true_divide(a, b, out=out, where=where)


def binary_search(function,
                  target,
                  lower_bound,
                  upper_bound,
                  tolerance=1e-4):
    lh = lower_bound
    rh = upper_bound
    while abs(rh - lh) > tolerance:
        mh = np.mean([lh, rh])
        lx, mx, rx = [function(h) for h in (lh, mh, rh)]
        if lx == target:
            return lx
        if rx == target:
            return rx

        if lx <= target and rx >= target:
            if mx > target:
                rh = mh
            else:
                lh = mh
        elif lx > target and rx < target:
            lh, rh = rh, lh
        else:
            return None
    return mh


def printProgressBar (iteration, total, prefix = '', suffix = '', decimals = 1, length = 100, fill = '█', printEnd = "\r"):
    """
    From shorturl.at/kpSY1
    Permanent link: shorturl.at/ADMU3

    Call in a loop to create terminal progress bar
    @params:
        iteration   - Required  : current iteration (Int)
        total       - Required  : total iterations (Int)
        prefix      - Optional  : prefix string (Str)
        suffix      - Optional  : suffix string (Str)
        decimals    - Optional  : positive number of decimals in percent complete (Int)
        length      - Optional  : character length of bar (Int)
        fill        - Optional  : bar fill character (Str)
        printEnd    - Optional  : end character (e.g. "\r", "\r\n") (Str)
    """
    percent = ("{0:." + str(decimals) + "f}").format(100 * (iteration / float(total)))
    filledLength = int(length * iteration // total)
    bar = fill * filledLength + '-' * (length - filledLength)
    print(f'\r{prefix} |{bar}| {percent}% {suffix}', end=printEnd)
    # Print New Line on Complete
    if iteration == total: 
        print()
