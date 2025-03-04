# ####################################################################### #
#                                                                         #
#   EEEEEE  XX  XX   AAAA   MM   MM  PPPPPP  LL      EEEEEE   11  55555   #
#   EE       XXXX   AA  AA  MMM MMM  PP  PP  LL      EE       11  55      #
#   EEEEE     XX    AA  AA  MMMMMMM  PPPPPP  LL      EEEEE    11  55555   #
#   EE       XXXX   AAAAAA  MM   MM  PP      LL      EE       11     55   #
#   EEEEEE  XX  XX  AA  AA  MM   MM  PP      LLLLLL  EEEEEE   11  55555   #
#                                                                         #
# ####################################################################### #

# Multiple criteria BMA model training: sea surface pressure ensemble
#  Vrugt, J.A. (2024), Distribution-Based Model Evaluation and            
#      Diagnostics: Elicitability, Propriety, and Scoring Rules for       
#      Hydrograph Functionals, Water Resources Research, 60,              
#      e2023WR036710, https://doi.org/10.1029/2023WR036710                
#  Vrugt, J.A., and B.A. Robinson (2007), Treatment of uncertainty using 
#      ensemble methods: Comparison of sequential data assimilation and 
#      Bayesian model averaging, Water Resources Research, 43, W01411, 
#      https://doi.org/10.1029/2005WR004838
#  Vrugt, J.A., M.P. Clark, C.G.H. Diks, Q. Duan, and B. A. Robinson      
#      (2006), Multi-objective calibration of forecast ensembles using     
#      Bayesian model averaging, Geophysical Research Letters, 33,        
#      L19817

import numpy as np
import os, sys
import pandas as pd
from scipy import io

current_dir = os.getcwd()                                       # Get the current working directory
parent_dir = os.path.abspath(os.path.join(current_dir, '..'))   # Go up one directory
from AMALGAM_BMA import setup_BMA
sys.path.append(parent_dir)                                     # add this to path
from AMALGAM import AMALGAM

# Define AMALGAM parameters
AMALGAMPar = {}
AMALGAMPar['N'] = 100  # Population size
AMALGAMPar['T'] = 250  # Number of generations
AMALGAMPar['m'] = 3    # Number of objective functions

# Define parameter information for the BMA model
Par_info = {}
Par_info['initial'] = 'latin'           # Latin hypercube sampling
Par_info['boundhandling'] = 'reflect'   # Boundary handling (reflective)

Func_name = 'AMALGAM_BMA.AMALGAM_BMA'   # Function name (model)

# Load the data
data = np.loadtxt('pressure.txt')

# Find the indices for the training period
id_1 = np.where((data[:, 0] == 2000) & (data[:, 1] == 4) & (data[:, 2] == 16))[0]
start_id = id_1[0]

id_2 = np.where((data[:, 0] == 2000) & (data[:, 1] == 6) & (data[:, 2] == 9))[0]
end_idx = id_2[-1]

# Define the training period (start and end indices)
T_idx = np.arange(start_id, end_idx + 1)  # Indices of the training period

# Extract ensemble forecasts and verifying observations
D = data[T_idx, 4:9]  # Ensemble forecasts (columns 5-9 in MATLAB, 4-8 in Python)
y = data[T_idx, 3]    # Verifying observations (column 4 in MATLAB, 3 in Python)

# Set up options for the BMA method
options = {}
options['BMA'] = 'yes'  # Activate BMA method
PDF = 'normal'          # Forecast PDF ('normal'/'gamma')
VAR = '2'               # Variance option ('1'/'2'/'3'/'4')

# Call the setup_BMA function to prepare the BMA model with bias correction
AMALGAMPar, Par_info, D_bc, A, B = setup_BMA(AMALGAMPar, Par_info, D, y, VAR)

# Create the plugin structure for BMA information
plugin = {}
plugin['BMA'] = {
    'PDF': PDF,
    'VAR': VAR,
    'D': D_bc,
    'y': y,
    'K': D.shape[1]  # Number of ensemble members (columns)
}

# Set options for the optimization process
options['print'] = 'yes'      # Print output to screen (figures)
options['save'] = 'yes'       # Save memory during AMALGAM restart run
options['modout'] = 'yes'     # Return simulation of BMA model

if __name__ == '__main__':
    # Run the AMALGAM code and obtain non-dominated solution set
    X, F, output, Z, YX = AMALGAM(AMALGAMPar, Func_name, Par_info, options, plugin)
