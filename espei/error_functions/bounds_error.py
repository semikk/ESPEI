"""
Defines an error that enforces bounds on the parameter values.
"""

import numpy as np

def generate_bounds(parameters, factor):
    """
    Generate (min, max) bounds for the given parameters.

    Parameters
    ----------
    parameters : list
        List of candidate parameter values
    factor : float or [float]
        Scalor or list of factor(s) to use in creating a min/max bound by p-f*p or p+f*p. If a list, it must be the same length as parameters.

    Returns
    -------
    list
        List of (min, max) tuples corresponding to each parameter

    Notes
    -----
    Bounds are generated for each parameter with the bounds being (x-f*x, x+f*x), where f is a scaling factor given by
    the factor argument. If factor is given as just a single scaling factor, it will be applied to all of the parameters.
    Factors can also be given as a list of scalars, corresponding to the scale factor, f, for each parameter individually.

    """

    if np.isscalar(factor):
        bounds = [(p-p*factor, p+p*factor) for p in parameters]
    else:
        bounds = [(p-p*f, p+p*f) for p, f in zip(parameters, factor)]
    return bounds


def calculate_bounds_error(parameters, bounds):
    """
    Return 0 (no error) if all parameters are inside the bounds and -inf otherwise.

    Parameters
    ----------
    parameters : list
        List of candidate parameter values
    bounds : list of tuples
        List of bound tuples corresponding to parameter (min, max). Must be the same length as parameters.

    Returns
    -------

    Notes

    This is different from a normalization error because it's a binary 0 or -inf. It's not weighted by parameter values.

    """
    bounds_violated = [(p < min_bound) or (p > max_bound) for p, (min_bound, max_bound) in zip(parameters, bounds)]
    if any(bounds_violated):
        return -np.inf
    else:
        return 0
