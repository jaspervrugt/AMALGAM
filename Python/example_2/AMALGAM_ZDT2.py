import numpy as np

def AMALGAM_ZDT2(x, d):
    """
    AMALGAM_ZDT2: ZDT2 bivariate test function of Zitzler et al. 2000
    
    Parameters:
    x : array-like
        1D parameter vector.
    d : int
        Number of parameters.

    Returns:
    Fx : ndarray
        Pair of objective function values.
    """
    
    f = x[0]                                # f(x)
    g = 1 + 9 / (d - 1) * np.sum(x[1:d])    # g(x)
    h = 1 - (f / g) ** 2                    # h(x)
    
    Fx = np.array([f, g * h])               # Pair of objective function values
    return Fx

