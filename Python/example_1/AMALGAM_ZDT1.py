import numpy as np

def AMALGAM_ZDT1(x, d):
    """
    AMALGAM_ZDT1: ZDT1 bivariate test function from Zitzler et al. 2000.

    Parameters:
    - x (numpy array): 1xd parameter vector
    - d (int): Number of parameters

    Returns:
    - Fx (numpy array): Pair of objective function values
    """
    f = x[0]                            # f(x)
    g = 1 + 9/(d - 1) * np.sum(x[1:d])  # g(x)
    h = 1 - np.sqrt(f / g)              # h(x)
    Fx = np.array([f, g * h])           # Pair of FX values
    z = np.random.uniform(0,1,(3))

    return Fx.flatten(), z.flatten()
