# ####################################################################### #
#                                                                         #
#   EEEEEE  XX  XX   AAAA   MM   MM  PPPPPP  LL      EEEEEE   11  222222  #
#   EE       XXXX   AA  AA  MMM MMM  PP  PP  LL      EE       11  22 22   #
#   EEEEE     XX    AA  AA  MMMMMMM  PPPPPP  LL      EEEEE    11    22    #
#   EE       XXXX   AAAAAA  MM   MM  PP      LL      EE       11   22     #
#   EEEEEE  XX  XX  AA  AA  MM   MM  PP      LLLLLL  EEEEEE   11  222222  #
#                                                                         #
# ####################################################################### #

# Multiple criteria hymod model training: driven & nondriven hydrograph
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

sys.path.append(os.path.join(parent_dir, 'miscellaneous'))	    # Add miscellaneous directory to Python path
from AMALGAM_functions import Compute_FX_true		            # Import functions

# Define AMALGAM parameters
AMALGAMPar = {
    'N': 100,  # Population size
    'T': 150,  # Generations
    'd': 5,    # Parameters
    'm': 2     # Objective functions
}

# Parameter Info
Par_info = {
    'initial': 'latin',                                 # Latin hypercube sampling
    'boundhandling': 'reflect',                         # Boundary handling
    'min': np.array([1.0, 0.10, 0.10, 0.00, 0.10]),     # Minimum parameter values
    'max': np.array([500, 2.00, 0.99, 0.10, 0.99]),     # Maximum parameter values
    'steps': np.full(AMALGAMPar['d'], 499),             # Steps for each parameter
    'names': ['C_max','b_exp','alfa','R_s','R_q']	
}

Func_name = 'AMALGAM_hymod.AMALGAM_hymod'           # Define name of function [explicit solution: not recommended]

# Load the data
bound = np.loadtxt('bound.txt')

# Constants
Area = 1944  # Area of Leaf River in km^2
conv_mult = Area * (1000 * 1000) / (1000 * 60 * 60 * 24)  # Convert m^3/s to mm/day
T_max = 795  # Use two years of data

# Extract data
Y_obs = bound[64:T_max, 3] / conv_mult  # Measured discharge in mm/day
PET = bound[0:T_max, 4]  # Potential Evapotranspiration (mm/day)
R = np.sum(bound[0:T_max, 5:9], axis=1)  # Daily rainfall (mm/day) (= sum of 6-hourly data)

# Indexes for driven and non-driven parts of the hydrograph
idx_d = np.where(R[64:] > 0)[0]  # Indexes for driven part
N_d = len(idx_d)  # Number of driven entries
idx_nd = np.where(R[64:] == 0)[0]  # Indexes for non-driven part
N_nd = len(idx_nd)  # Number of non-driven entries

# Plugin structure (used to pass to AMALGAM)
plugin = {
    'fieldnames': ['T_max', 'Y_obs', 'PET', 'R', 'idx_d', 'N_d', 'idx_nd', 'N_nd'],
    'T_max': T_max,
    'Y_obs': Y_obs,
    'PET': PET,
    'R': R,
    'idx_d': idx_d,
    'N_d': N_d,
    'idx_nd': idx_nd,
    'N_nd': N_nd
}

# Options for AMALGAM
options = {
    'print': 'yes',  	# Print output to screen (figures)
    'parallel': 'no',  	# Multi-core: No IO writing
    'IO': 'no',  	# No I/O writing workers (default)
    'modout': 'yes',  	# Return hmodel simulations
}

if __name__ == '__main__':
    # Run the AMALGAM code and obtain non-dominated solution set
    X, F, output, Z, YX = AMALGAM(AMALGAMPar, Func_name, Par_info, options, plugin)
