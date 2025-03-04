import numpy as np

def AMALGAM_ZDT4(x, d):
    """
    AMALGAM_ZDT4: ZDT4 bivariate test function of Zitzler et al. 2000.
    
    Parameters:
    - x (ndarray): 1d parameter vector (length d)
    - d (int): Number of parameters
    
    Returns:
    - Fx (ndarray): Pair of objective function values [f(x), g(x)*h(x)]
    """
    f = x[0]                                                                    # f(x) is simply the first element of x
    g = 1 + 10 * (d - 1) + np.sum(x[1:]**2 - 10 * np.cos(4 * np.pi * x[1:]))    # g(x)
    h = 1 - np.sqrt(f / g)                                                      # h(x)
    
    Fx = np.array([f, g * h])                                                   # Pair of objective function values
    
    return Fx
