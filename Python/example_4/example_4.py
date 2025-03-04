# ####################################################################### #
#                                                                         #
#   EEEEEE  XX  XX   AAAA   MM   MM  PPPPPP  LL      EEEEEE        4444   #
#   EE       XXXX   AA  AA  MMM MMM  PP  PP  LL      EE           44 44   #
#   EEEEE     XX    AA  AA  MMMMMMM  PPPPPP  LL      EEEEE       44  44   #
#   EE       XXXX   AAAAAA  MM   MM  PP      LL      EE          444444   #
#   EEEEEE  XX  XX  AA  AA  MM   MM  PP      LLLLLL  EEEEEE          44   #
#                                                                         #
# ####################################################################### #

# ZDT3 example from the following paper
#  Zitzler, E., K. Deb, and L. Thiele (2000), Comparison of Multiobjective
#      Evolutionary Algorithms: Empirical Results, Evolutionary 
#      Computation, 8 (2), 183-195, 2000
#      https://sop.tik.ee.ethz.ch/publicationListFiles/zdt2000a.pdf

import numpy as np
import os, sys

current_dir = os.getcwd()                                       # Get the current working directory
parent_dir = os.path.abspath(os.path.join(current_dir, '..'))   # Go up one directory
sys.path.append(parent_dir)                                     # add this to path
from AMALGAM import AMALGAM

sys.path.append(os.path.join(parent_dir, 'miscellaneous'))	    # Add miscellaneous directory to Python path
from AMALGAM_functions import Compute_FX_true		            # Import functions

# Define AMALGAM parameters
AMALGAMPar = {}
AMALGAMPar['N'] = 100    # Population size
AMALGAMPar['T'] = 100    # Number of generations
AMALGAMPar['d'] = 10     # Number of parameters
AMALGAMPar['m'] = 2      # Number of objective functions

# Parameter information for Latin hypercube sampling-based initialization
Par_info = {}
Par_info['initial'] = 'latin'                                               # Latin hypercube sampling
Par_info['boundhandling'] = 'bound'                                         # Boundary handling
Par_info['min'] = np.concatenate(([0], -5 * np.ones(AMALGAMPar['d'] - 1)))  # Minimum parameter values
Par_info['max'] = np.concatenate(([1], 5 * np.ones(AMALGAMPar['d'] - 1)))   # Maximum parameter values

# Define the function name
Func_name = 'AMALGAM_ZDT4.AMALGAM_ZDT4'  # The benchmark function

# Plugin for the second argument of Func_name (number of parameters)
plugin = AMALGAMPar['d']

# Set optimization options
options = {}
options['print'] = 'yes'   # Print output to screen (figures)
options['parallel'] = 'no'

# Compute the true Pareto front for the DTLZ2 problem
FX_true = Compute_FX_true(AMALGAMPar, Func_name, 'ZDT4', 500, plugin)

if __name__ == '__main__':
    # Run the AMALGAM optimization and obtain non-dominated solution set
    X, FX, output, Z, YX = AMALGAM(AMALGAMPar, Func_name, Par_info, options, plugin, FX_true)
