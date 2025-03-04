import numpy as np

def AMALGAM_DTLZ1(x, M):
    """
    AMALGAM_DTLZ1: DTLZ1 m-variate test function of Deb et al. 2001.
    
    Parameters:
    - x (ndarray): 1d parameter vector (length d)
    - M (int): Number of objectives
    
    Returns:
    - Fx (ndarray): m objective function values
    """
    M = int(M)
    d = len(x)      # Number of decision variables
    k = d - M + 1   # Integer k

    # Calculate g (the "decoupled" part of the DTLZ1 function)
    g = 100 * (k + np.sum((x[M-1:] - 0.5)**2 - np.cos(20 * np.pi * (x[M-1:] - 0.5))))
    # Initialize the objective vector Fx
    Fx = np.full(M, np.nan)
    # Compute the first objective
    Fx[0] = 0.5 * np.prod(x[:M-1]) * (1 + g)
    # Compute the remaining objectives
    for m in range(1, M-1):
        Fx[m] = 0.5 * np.prod(x[:M-m-1]) * (1 - x[M-m-1]) * (1 + g)

    # Compute the last objective
    Fx[M-1] = 0.5 * (1 - x[0]) * (1 + g)
    
    return Fx

