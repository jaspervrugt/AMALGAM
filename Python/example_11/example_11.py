# ####################################################################### #
#                                                                         #
#   EEEEEE  XX  XX   AAAA   MM   MM  PPPPPP  LL      EEEEEE     11  11    #
#   EE       XXXX   AA  AA  MMM MMM  PP  PP  LL      EE         11  11    #
#   EEEEE     XX    AA  AA  MMMMMMM  PPPPPP  LL      EEEEE      11  11    #
#   EE       XXXX   AAAAAA  MM   MM  PP      LL      EE         11  11    #
#   EEEEEE  XX  XX  AA  AA  MM   MM  PP      LLLLLL  EEEEEE     11  11    #
#                                                                         #
# ####################################################################### #

# Hydrologic modeling using hmodel: Driven & nondriven hydrograph
#  Vrugt, J.A., H.V. Gupta, L.A. Bastidas, W. Bouten, and S. Sorooshian
#      (2000), Effective and efficient algorithm for multiobjective 
#      optimization of hydrologic models, Water Resources Research, 
#      39 (8), 1214, https://doi.org./10.1029/2002WR001746

import numpy as np
import os, sys
import pandas as pd
from scipy import io

current_dir = os.getcwd()                                       # Get the current working directory
parent_dir = os.path.abspath(os.path.join(current_dir, '..'))   # Go up one directory
sys.path.append(parent_dir)                                     # add this to path
from AMALGAM import AMALGAM

# Define AMALGAM parameters
AMALGAMPar = {
    'N': 50,   # Population size
    'T': 50,   # Generations
    'd': 7,    # Parameters
    'm': 2     # Objective functions
}

# Parameter Info
Par_info = {
    'initial': 'latin',                                         # Latin hypercube sampling
    'boundhandling': 'reflect',                                 # Boundary handling
    'min': np.array([1, 10, 0.1, 0.1, -10, 0.1, 0.1]),          # Minimum parameter values
    'max': np.array([10, 1000, 100, 100, 10, 10, 150]),         # Maximum parameter values
    'names': ['I_max', 'S_u,max', 'Q_s,max', 'alpha_E', 'alpha_F', 'K_f', 'K_s']
}

Func_name = 'AMALGAM_hmodel.AMALGAM_hmodel'                     # Define name of function

# Data multi-criteria hmodel training 
mopex = np.loadtxt('03451500.dly')

# Define training data set
idx = np.where((mopex[:, 0] > 1959) & (mopex[:, 0] < 1999))[0]
n = mopex.shape[0]
tout = np.arange(n+1)           # Time vector for observations

Y_obs = mopex[idx[0:n], 5]      # Measured discharge data (mm/day)
Y_obs = Y_obs[730:n]            # Adjust the length as in the original code

# Data variables
data = {
    'P': mopex[idx[0:n], 3],    # Daily rainfall (mm/d)
    'Ep': mopex[idx[0:n], 4],   # Daily evaporation (mm/d)
    'aS': 1e-6                  # Percolation coefficient
}

# hmodel options
hmodel_opt = {
    'InitialStep': 1,   # Initial time-step (d)
    'MaxStep': 1,       # Maximum time-step (d)
    'MinStep': 1e-6,    # Minimum time-step (d)
    'RelTol': 1e-3,     # Relative tolerance
    'AbsTol': 1e-3,     # Absolute tolerances (mm)
    'Order': 2          # 2nd order method (Heun)
}

y0 = 1e-6 * np.ones(5)  # Initial conditions

# Index for driven and non-driven parts
id_d = np.where(data['P'][730:n] > 0)[0]
N_d = len(id_d)         # observations driven part    

id_nd = np.where(data['P'][730:n] == 0)[0]
N_nd = len(id_nd)       # observations non-driven part

# Plugin structure (used to pass to the AMALGAM algorithm)
plugin = {'fields': ['fieldNames', 'tout', 'data', 'hmodel_opt', 'y0', 'Y_obs', 'n', 'id_d', 'N_d', 'id_nd', 'N_nd'],
          'fieldNames': ['fieldNames', 'tout', 'data', 'hmodel_opt', 'y0', 'Y_obs', 'n', 'id_d', 'N_d', 'id_nd', 'N_nd'],
          'tout': tout,
          'data': data,
          'hmodel_opt': hmodel_opt,
          'y0': y0,
          'Y_obs': Y_obs,
          'n': n,
          'id_d': id_d,
          'N_d': N_d,
          'id_nd': id_nd,
          'N_nd': N_nd}

# Options for AMALGAM
options = {
    'print': 'yes',      # Print output to screen (figures)
    'parallel': 'yes',   # Multi-core: No IO writing
    'modout': 'yes',     # Return hmodel simulations
    'save': 'yes'        # Save memory during trial
}

if __name__ == '__main__':
    # Run the AMALGAM code and obtain non-dominated solution set
    X, FX, output, Z, YX = AMALGAM(AMALGAMPar, Func_name, Par_info, options, plugin)
