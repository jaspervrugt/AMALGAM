import numpy as np

############## HYMOD: Explicit Euler [= not recommended]
def AMALGAM_hymod(par, plugin):
    # Unpack parameters
    cmax, bexp, alfa, Rs, Rq = par
    T_max, Y_obs, PET, R, idx_d, N_d, idx_nd, N_nd = plugin['T_max'], plugin['Y_obs'], plugin['PET'], plugin['R'], plugin['idx_d'], plugin['N_d'], plugin['idx_nd'], plugin['N_nd']
    
    # Initialize states
    x_l = 0.0  # Initial state loss tank
    x_s = 0     # Initial state slow tank
    x_q = np.zeros(3)  # Initial states fast tanks
    output = np.full(T_max, np.nan)  # Simulated discharge
    
    # Main loop for simulation
    for t in range(T_max):
        # Compute excess precipitation and evaporation
        ER1, ER2, x_l = excess(x_l, cmax, bexp, R[t], PET[t])
        # Calculate total effective rainfall
        ET = ER1 + ER2
        # Partition ER between quick and slow flow reservoirs
        UQ = alfa * ET
        US = (1 - alfa) * ET
        # Route slow flow component with a single linear reservoir
        x_s, QS = linres(x_s, US, Rs)
        # Route quick flow component with linear reservoirs
        q_in = UQ
        for k in range(3):
            x_q[k], q_out = linres(x_q[k], q_in, Rq)
            q_in = q_out
        # Compute total flow for timestep (in mm/day)
        output[t] = QS + q_out
    
    # Apply burn-in of 65 days
    Y_sim = output[64:T_max]
    # Compute RMSE for driven and non-driven parts
    F = np.zeros(2)
    F[0] = np.sqrt(np.sum((Y_sim[idx_d] - Y_obs[idx_d])**2) / N_d)  # RMSE driven part
    F[1] = np.sqrt(np.sum((Y_sim[idx_nd] - Y_obs[idx_nd])**2) / N_nd)  # RMSE nondriven part
    
    return F, Y_sim


# Secondary functions

def excess(x_loss, cmax, bexp, Pval, PETval):
    xn_prev = x_loss
    ct_prev = cmax * (1 - ((1 - ((bexp + 1) * xn_prev / cmax)) ** (1 / (bexp + 1))))
    # Calculate effective rainfall 1
    ER1 = max((Pval - cmax + ct_prev), 0.0)
    Pval -= ER1
    dummy = min(((ct_prev + Pval) / cmax), 1)
    xn = (cmax / (bexp + 1)) * (1 - (1 - dummy)**(bexp + 1))
    # Calculate effective rainfall 2
    ER2 = max(Pval - (xn - xn_prev), 0)
    # Alternative approach for evapotranspiration
    evap = (1 - ((cmax / (bexp + 1)) - xn) / (cmax / (bexp + 1))) * PETval
    xn = max(xn - evap, 0)
    
    return ER1, ER2, xn


def linres(x, q_in, R):
    # Linear reservoir function
    x = (1 - R) * x + (1 - R) * q_in
    q_out = (R / (1 - R)) * x
    return x, q_out
