import numpy as np

def AMALGAM_ZDT3(x, d):
    """
    AMALGAM_ZDT3: ZDT3 bivariate test function of Zitzler et al. 2000.
    
    Parameters:
    - x (ndarray): 1d parameter vector (length d)
    - d (int): Number of parameters
    
    Returns:
    - Fx (ndarray): Pair of objective function values [f(x), g(x)*h(x)]
    """
    f = x[0]                                                    # f(x) is simply the first element of x
    g = 1 + (9 / (d - 1)) * np.sum(x[1:])                       # g(x) is the sum of the remaining elements
    h = 1 - np.sqrt(f / g) - (f / g) * np.sin(10 * np.pi * f)   # h(x) term
    
    Fx = np.array([f, g * h])                                   # Pair of objective function values
    
    return Fx
