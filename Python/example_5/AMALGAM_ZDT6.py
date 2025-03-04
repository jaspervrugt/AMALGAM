import numpy as np

def AMALGAM_ZDT6(x, d):
    """
    AMALGAM_ZDT6: ZDT6 bivariate test function of Zitzler et al. 2000.
    
    Parameters:
    - x (ndarray): 1d parameter vector (length d)
    - d (int): Number of parameters
    
    Returns:
    - Fx (ndarray): Pair of objective function values [f(x), g(x)*h(x)]
    """
    f = 1 - np.exp(-4 * x[0]) * np.sin(6 * np.pi * x[0])**6                     # f(x)
    g = 1 + 9 * (1 / (d - 1) * np.sum(x[1:]))**(1 / 4)                          # g(x)
    h = 1 - (f / g)**2                                                          # h(x)
    
    Fx = np.array([f, g * h])                                                   # Pair of objective function values
    
    return Fx
