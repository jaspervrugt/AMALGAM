import numpy as np

def AMALGAM_DTLZ3(x, M):
    """
    AMALGAM_DTLZ3: DTLZ3 m-variate test function of Deb et al. 2001.
    
    Parameters:
    - x (ndarray): 1D parameter vector (length d)
    - M (int): Number of objectives
    
    Returns:
    - Fx (ndarray): m objective function values
    """
    Fx = np.full(M, np.nan)     # Initialize the objective values array
    d = len(x)                  # Number of decision variables
    k = d - M + 1               # Integer k
    xM = x[-k:]                 # Last k variables: xM
    xnotM = x[:-k]              # Remaining parameters
    
    #    k1 = d - k  # Number of remaining decision variables
    k1 = d-k    

    # Calculate g(xM)
    g = 100 * (k + np.sum((xM - 0.5)**2 - np.cos(20 * np.pi * (xM - 0.5))))
    
    # Calculate the first objective f_1(x)
    Fx[0] = (1 + g) * np.prod(np.cos(xnotM[:k1] * np.pi / 2))
    
    # Calculate the intermediate objectives f_2(x) to f_{M-1}(x)
    for m in range(1, M - 1):
        Fx[m] = (1 + g) * np.prod(np.cos(xnotM[:k1 - m] * np.pi / 2)) * np.sin(xnotM[k1 - m] * np.pi / 2)
    
    # Calculate the last objective f_M(x)
    Fx[M - 1] = (1 + g) * np.sin(xnotM[0] * np.pi / 2)
    
    return Fx
