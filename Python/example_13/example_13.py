# ####################################################################### #
#                                                                         #
#   EEEEEE  XX  XX   AAAA   MM   MM  PPPPPP  LL      EEEEEE   11  33333   #
#   EE       XXXX   AA  AA  MMM MMM  PP  PP  LL      EE       11     33   #
#   EEEEE     XX    AA  AA  MMMMMMM  PPPPPP  LL      EEEEE    11    333   #
#   EE       XXXX   AAAAAA  MM   MM  PP      LL      EE       11     33   #
#   EEEEEE  XX  XX  AA  AA  MM   MM  PP      LLLLLL  EEEEEE   11  33333   #
#                                                                         #
# ####################################################################### #

# Multiple criteria BMA model training: river discharge ensemble
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

# Set up the parameters for AMALGAM
AMALGAMPar = {}
AMALGAMPar['N'] = 100  # Population size
AMALGAMPar['T'] = 100  # Number of generations
AMALGAMPar['m'] = 3    # Number of objective functions

# Set up parameter info for the BMA model
Par_info = {}
Par_info['initial'] = 'latin'           # Latin hypercube sampling
Par_info['boundhandling'] = 'reflect'   # Explicit boundary handling

Func_name = 'AMALGAM_BMA.AMALGAM_BMA'   # Define the name of the function

# Load the data
data = np.loadtxt('discharge.txt')

# Define the training period and select the corresponding data
T_idx = np.arange(0, 3000)              # Indices for the training period (0-based)
D = data[T_idx, :8]                     # Ensemble forecasts (assuming 8 models in columns)
y = data[T_idx, 8]                      # Verifying observations (9th column in MATLAB, 8th in Python)

# Set up options for BMA method
options = {}
options['BMA'] = 'yes'  # Activate BMA method
PDF = 'normal'          # Forecast pdf ('normal'/'gamma')
VAR = '4'               # Variance option ('1'/'2'/'3'/'4')

# Call the setup_BMA function to prepare the BMA model with bias correction
AMALGAMPar, Par_info, D_bc, A, B = setup_BMA(AMALGAMPar, Par_info, D, y, VAR)

# Structure for BMA info
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
options['modout'] = 'yes'     # Return hmodel simulations
options['save'] = 'yes'       # Save memory during AMALGAM restart run

if __name__ == '__main__':
    # Run the AMALGAM code (you would need to have a Python equivalent for AMALGAM)
    X, F, output, Z, YX = AMALGAM(AMALGAMPar, Func_name, Par_info, options, plugin)
