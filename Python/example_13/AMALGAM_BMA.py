import numpy as np
from scipy.stats import (norm, gamma)

## BMA case study
def AMALGAM_BMA(par, plugin):
    """
    This function calculates the RMSE, Ignorance Score (IS), and Continuous Ranked Probability Score (CRPS)
    of the Bayesian Model Averaging (BMA) mixture distribution.
    
    Parameters:
        par (array): Parameters, including the weights for the ensemble members and variances.
        plugin (dict): A dictionary containing the necessary data, such as ensemble forecasts and verifying data.
    
    Returns:
        FX (array): Array containing the RMSE, Ignorance Score, and CRPS.
        G_dot (array): The BMA forecast mean.
    """
    
    D = plugin['BMA']['D']  # Ensemble forecasts
    y = plugin['BMA']['y']  # Verifying data
    n, K = D.shape          # Number of forecasts (n) and number of ensemble members (K)
    
    # Unpack weights
    w = par[:K]
    
    # Variance options
    VAR = plugin['BMA']['VAR']
    
    if VAR == '1':      # Common constant variance
        sigma = par[K] * np.ones((n, K))
    elif VAR == '2':    # Individual constant variance
        sigma = np.multiply(par[K:2*K], np.ones((n, K)))
    elif VAR == '3':    # Common non-constant variance
        c = par[K]
        sigma = c * D
    elif VAR == '4':    # Individual non-constant variance
        c = par[K:2*K]
        sigma = np.multiply(c, D)
    
    sigma = np.maximum(sigma, np.finfo(float).eps)  # Ensure sigma is not zero
    
    # Conditional distribution options
    PDF = plugin['BMA']['PDF']
    
    # Calculate likelihoods
    Y = np.tile(y, (K, 1)).T  # Make K copies of verifying data
    
    if PDF == 'normal':  # Gaussian distribution
        A = D
        B = sigma
        L = norm.pdf(Y, A, B)

    elif PDF == 'gamma':  # Gamma distribution
        mu = np.abs(D)
        var = sigma ** 2
        A = mu ** 2 / var
        B = var / mu
        L = gamma.pdf(Y, A, scale=B)

    # Calculate likelihoods
#    Y = np.tile(y, (K, 1)).T  # Make K copies of verifying data
#    L = pdf(PDF, Y, A, B)  # Assuming pdf function is defined elsewhere
    lik = np.sum(L * w, axis=1) + np.finfo(float).tiny  # BMA likelihoods
    
    G_dot = np.dot(D, w)  # BMA mean forecast
    res = y - G_dot  # Residuals
    
    # Generate two samples for CRPS calculation
    G = BMA_rnd(PDF, w, A, B, 2)  # Assuming BMA_rnd function is defined elsewhere
    
    # Calculate RMSE, Ignorance Score, and CRPS
    RMSE = np.sqrt(np.sum(res ** 2) / n)  # RMSE of the mean forecast
    IS = -np.mean(np.log(lik))  # Ignorance Score
    CRPS = np.mean(np.abs(G[0] - y)) - 0.5 * np.mean(np.abs(G[0] - G[1]))  # CRPS score - negatively oriented
    
    FX = np.array([RMSE, IS, CRPS])  # Return RMSE, IS, and CRPS as an array
    
    return FX, G_dot


def setup_BMA(AMALGAMPar, Par_info, D, y, VAR):
    """
    Setup BMA mixture model estimation
    
    Parameters:
    - AMALGAMPar: Dictionary containing parameters
    - Par_info: Dictionary to store parameter information
    - D: NumPy array of forecasts (shape: n x K)
    - y: NumPy array of verifying data (shape: n x 1)
    - VAR: String indicating variance type ('1', '2', '3', or '4')

    Returns:
    - AMALGAMPar: Updated dictionary with new parameters
    - Par_info: Updated dictionary with parameter information
    - D_bc: Bias-corrected forecasts
    - A: Intercepts for bias correction
    - B: Slopes for bias correction
    """
    
    # Get number of forecasts (n) and number of ensemble members (K)
    n, K = D.shape
    
    adjust = True  # Linear bias correction
    
    # Perform bias correction using the ComputeAB function
    D_bc, A, B = ComputeAB(D, y, n, K, adjust)
    
    # Initialize parameter names
    par_name = [f'\\beta_{z+1}' for z in range(K)]  # Beta names for the K ensemble members
    
    # Set parameter names and adjust settings based on VAR value
    if VAR == '1':  # Common constant variance
        AMALGAMPar['d'] = K + 1
        Par_info['max'] = np.hstack([np.ones(K), 2 * np.std(y)])
        par_name.append('\\sigma')
    elif VAR == '2':  # Individual constant variance
        AMALGAMPar['d'] = 2 * K
        Par_info['max'] = np.hstack([np.ones(K), 2 * np.std(y) * np.ones(K)])
        for z in range(K):
            par_name.append(f'\\sigma_{z+1}')
    elif VAR == '3':  # Common non-constant variance
        AMALGAMPar['d'] = K + 1
        Par_info['max'] = np.hstack([np.ones(K), 2])
        par_name.append('c')
    elif VAR == '4':  # Individual non-constant variance
        AMALGAMPar['d'] = 2 * K
        Par_info['max'] = np.hstack([np.ones(K), 2 * np.ones(K)])
        for z in range(K):
            par_name.append(f'c_{z+1}')
    else:
        raise ValueError("Unknown variance option; choose between '1', '2', '3', or '4'")

    # Final adjustments to Par_info
    Par_info['names'] = par_name  # Set parameter names
    Par_info['unit_simplex'] = 'yes'  # Weights on unit Simplex
    Par_info['min'] = np.zeros(AMALGAMPar['d'])  # Min. values for BMA weights/vars
    
    return AMALGAMPar, Par_info, D_bc, A, B

# ComputeAB function (helper for bias correction)
def ComputeAB(D, y, n, K, adjust=True):
    """
    Performs linear bias correction if indicated
    
    Parameters:
    - D: NumPy array of ensemble forecasts (shape: n x K)
    - y: NumPy array of verifying data (shape: n x 1)
    - n: Number of forecasts
    - K: Number of ensemble members
    - adjust: Boolean flag to indicate if bias correction should be applied

    Returns:
    - D_bc: Bias-corrected ensemble forecasts
    - A: Intercepts of the linear bias correction
    - B: Slopes of the linear bias correction
    """
    
    if adjust:
        B = np.zeros(K)
        A = np.zeros(K)
        for k in range(K):
            T = np.cov(D[:, k], y)[0, 1]  # Covariance
            B[k] = T / np.var(D[:, k])  # Slope of the linear regression
            A[k] = np.mean(y) - B[k] * np.mean(D[:, k])  # Intercept of the regression
    else:
        B = np.ones(K)
        A = np.zeros(K)
    
    # Compute the bias-corrected forecasts
    D_bc = np.full_like(D, np.nan)
    for k in range(K):
        D_bc[:, k] = B[k] * D[:, k] + A[k]
    
    return D_bc, A, B


def BMA_rnd(PDF, w, A, B, P):
    """
    Samples trajectories from the BMA forecast distribution.

    Parameters:
        PDF (str): The type of the probability distribution to sample from ('normal' or 'gamma').
        w (array): The weights for each ensemble member (length K).
        A (array): The shape (n, K) matrix of distribution parameters (e.g., mean or shape).
        B (array): The shape (n, K) matrix of distribution parameters (e.g., standard deviation or scale).
        P (int): The number of trajectories to sample.

    Returns:
        G (list): A list containing P arrays of sampled trajectories from the BMA distribution.
    """
    n, K = A.shape  # Number of forecasts (n) and ensemble members (K)

    A = np.asarray(A)
    B = np.asarray(B)

    G = []  # List to hold the sampled trajectories

    for p in range(P):
        # Draw n samples from the ensemble members (1 to K) with weights w
        id = np.random.choice(K, n, p=w)  # Weighted sampling
        G_p = np.zeros(n)  # Array to hold the sampled trajectory for this p-th sample
        for k in range(K):
            # Find indices of samples corresponding to the k-th mixture component
            T = np.where(id == k)[0]
            # Sample from the k-th mixture component
            if PDF == 'normal':
                if len(T) > 0:  # Check if T is not empty
                    G_p[T] = np.random.normal(A[T, k], B[T, k])  # Normal distribution
            elif PDF == 'gamma':
                if len(T) > 0:  # Check if T is not empty
                    G_p[T] = np.random.gamma(A[T, k], B[T, k])  # Gamma distribution
            else:
                raise ValueError("Unsupported PDF type. Use 'normal' or 'gamma'.")

        G.append(G_p)
    
    return G


def pdf(name, x, *params):
    """
    Compute the PDF (Probability Density Function) for a specified distribution.
    
    Parameters:
        name (str): The name of the distribution (e.g., 'beta', 'normal', etc.).
        x (array-like): The values at which to evaluate the PDF.
        *params (list): The parameters required by the distribution.
        
    Returns:
        y (array-like): The values of the PDF evaluated at the points in x.
    """
    # Check the distribution name and compute the PDF accordingly
    if name in ['gamma', 'Gam']:
        a, scale = params
        y = gamma.pdf(x, a, scale=scale)
        
    elif name in ['normal', 'Norm']:
        mean, std = params
        y = norm.pdf(x, mean, std)
        
    else:
        raise ValueError(f"Unsupported distribution: {name}")
        
    return y
