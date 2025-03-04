import numpy as np

############### HMODEL ####################
def AMALGAM_hmodel(par, plugin):
    """
    hmodel simulation of discharge according to Schoups et al. 2010
    Returns Root Mean Square Error (RMSE) of driven and non-driven part hydrograph.
    
    Parameters:
        par (array): Parameter values for the model.
        plugin (dict): Dictionary containing plugin information, including time vector, model options, observed data, etc.
    
    Returns:
        F (array): RMSE of driven and non-driven part hydrograph.
        Y_sim (array): Simulated discharge after burn-in period.
    """
    
    # Extract various input data from the plugin dictionary
    tout = plugin['tout']
    data = plugin['data']
    hmodel_opt = plugin['hmodel_opt']
    y0 = plugin['y0']
    Y_obs = np.array(plugin['Y_obs']).astype(float)     # vector of reals
    n = int(plugin['n'])                                # scalar
    id_d = np.array(plugin['id_d']).astype(int)         # vector of integers
    N_d = int(plugin['N_d'])                            # scalar
    id_nd = np.array(plugin['id_nd']).astype(int)       # vector of integers
    N_nd = int(plugin['N_nd'])                          # scalar
    
    # Simulate discharge using the hmodel function
    y = hmodel(par, tout, data, hmodel_opt, y0)  # Assuming hmodel is defined elsewhere
    
    # Compute discharge from state
    #Y = y[4, 1:n] - y[4, 0:n-1]
    Y = y[4, 1:n+1] - y[4, 0:n]
 
     # Apply burn-in (skip first 731 data points, assumed to be 2 years)
    Y_sim = Y[730:n]

    # Calculate RMSE for driven part of hydrograph
    RMSE_driven = np.sqrt(np.sum((Y_sim[id_d] - Y_obs[id_d]) ** 2) / N_d)
    
    # Calculate RMSE for non-driven part of hydrograph
    RMSE_nondriven = np.sqrt(np.sum((Y_sim[id_nd] - Y_obs[id_nd]) ** 2) / N_nd)
    
    # Return the RMSE values and the simulated discharge
    F = np.array([RMSE_driven, RMSE_nondriven])

    return F, Y_sim


def hmodel(x, tout, data, options, y0):
    """
    Runs the hmodel and returns the driven and non-driven part.
    
    This Python version assumes that the crr_model function is implemented 
    (either as a C extension or as Python code). The provided parameters are:
    x      : vector of model parameters (e.g., interception, storage, etc.)
    tout   : time points for the integration
    data   : model data, dictionary with required fields
    options: options for the model (e.g., tolerance, step size, etc.)
    y0     : initial conditions for the model
    """
    
    # Update data dictionary with parameters from vector x
    data['Imax'] = x[0]     # interception storage capacity (mm)
    data['Sumax'] = x[1]    # unsaturated zone storage capacity (mm)
    data['Qsmax'] = x[2]    # maximum percolation rate (mm/d)
    data['aE'] = x[3]       # evaporation coefficient
    data['aF'] = x[4]       # runoff coefficient
    data['aS'] = 1e-6       # percolation coefficient (constant)
    data['Kf'] = x[5]       # fast-flow response time (d)
    data['Ks'] = x[6]       # slow-flow response time (d)
    
    # Run model function (assumed to be implemented in Python or using C extension)
    y = crr_model(tout, y0, data, options)
    
    return y


# Example of how this might be used
def crr_model(tout, y0, data, options):
    nvar = len(y0)  # Number of state variables
    nt = len(tout)  # Number of time steps
    
    # Call the Runge-Kutta solver
    return runge_kutta(nvar, nt, tout, y0, data, options)


def runge_kutta(nvar, nt, tout, y0, data, options):
    hin = options['InitialStep']
    hmax_ = options['MaxStep']
    hmin_ = options['MinStep']
    reltol = options['RelTol']
    abstol = options['AbsTol']
    order = options['Order']
    
    LTE = np.zeros(nvar)
    ytmp = np.zeros(nvar)
    w = np.zeros(nvar)
    ns = nt - 1

    # Initialize y
    y = np.zeros((nvar, nt))
    y[:, 0] = y0
    
    for s in range(1, ns + 1):
        t1 = tout[s-1]
        t2 = tout[s]
        
        h = hin
        h = max(hmin_, min(h, hmax_))
        h = min(h, t2 - t1)
        
        y[:, s] = y[:, s - 1]  # Set initial y
        
        t = t1
        while t < t2:
            ytmp[:] = y[:, s]  # Copy current y
            
            # RK2 integration step
            rk2(data, nvar, s, t, h, ytmp, LTE)
            
            # Check if the step is acceptable
            accept = 0
            wrms = 0
            for i in range(nvar):
                w[i] = 1.0 / (reltol * abs(ytmp[i]) + abstol)
                wrms += (w[i] * LTE[i])**2
            wrms = np.sqrt(wrms / nvar)
            if wrms <= 1:
                accept = 1
            
            if accept > 0:
                y[:, s] = ytmp  # Update y
                t += h
            
            # Adjust step size
            h = h * max(0.2, min(5.0, 0.9 * wrms**(-1.0/order)))
            h = max(hmin_, min(h, hmax_))
            h = min(h, t2 - t)
    
    return y


def rk2(data, nvar, s, t, h, u, LTE):
    udotE = np.zeros(nvar)
    uE = np.zeros(nvar)
    udot = np.zeros(nvar)

    # Euler solution
    flag = fRhs(s, t, u, udotE, data)
    uE[:] = u + h * udotE
    # Heun solution
    flag = fRhs(s, t + h, uE, udot, data)
    u[:] = u + 0.5 * h * (udotE + udot)
    
    # Estimate LTE (Local Truncation Error)
    LTE[:] = np.abs(uE - u)


def fRhs(s, t, u, udot, data):
    P = data['P']
    Ep = data['Ep']
    Imax = data['Imax']
    Sumax = data['Sumax']
    Qsmax = data['Qsmax']
    aE = data['aE']
    aF = data['aF']
    aS = data['aS']
    Kf = data['Kf']
    Ks = data['Ks']

    Si = u[0]
    Su = u[1]
    Sf = u[2]
    Ss = u[3]

    Precip = P[s-1]
    
    if Imax > 0.0:
        EvapI = Ep[s-1] * expFlux(Si / Imax, 50.0)
        P_e = P[s-1] * expFlux(Si / Imax, -50.0)
        Ep_e = max(0.0, Ep[s-1] - EvapI)
    else:
        EvapI = 0.0
        P_e = P[s-1]
        Ep_e = Ep[s-1]

    Evap = Ep_e * expFlux(Su / Sumax, aE)
    Perc = Qsmax * expFlux(Su / Sumax, aS)
    Runoff = P_e * expFlux(Su / Sumax, aF)
    FastQ = Sf / Kf
    SlowQ = Ss / Ks

    udot[0] = Precip - EvapI - P_e
    udot[1] = P_e - Evap - Perc - Runoff
    udot[2] = Runoff - FastQ
    udot[3] = Perc - SlowQ
    udot[4] = FastQ + SlowQ

    return 0


def expFlux(Sr, a):
    Sr = max(0.0, min(1.0, Sr))
    if abs(a) < 1e-6:
        Qr = Sr  # Linear approximation
    else:
        Qr = (1.0 - exponen(-a * Sr)) / (1.0 - exponen(-a))
    return Qr


def exponen(x):
    return np.exp(min(300.0, x))
