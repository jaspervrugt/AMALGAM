# ####################################################################### #
#                                                                         #
#   EEEEEE  XX  XX   AAAA   MM   MM  PPPPPP  LL      EEEEEE      888888   #
#   EE       XXXX   AA  AA  MMM MMM  PP  PP  LL      EE          88  88   #
#   EEEEE     XX    AA  AA  MMMMMMM  PPPPPP  LL      EEEEE       888888   #
#   EE       XXXX   AAAAAA  MM   MM  PP      LL      EE          88  88   #
#   EEEEEE  XX  XX  AA  AA  MM   MM  PP      LLLLLL  EEEEEE      888888   #
#                                                                         #
# ####################################################################### #

# DTLZ2 example from the following paper
#  Deb, K., L. Thiele, M. Laumanns, and E. Zitzler (2001), Scalable Test  
#      Problems for Evolutionary Multi-Objective Optimization. Kanpur,    
#      India: Kanpur Genetic Algorithms Lab. (KanGAL), Indian Institute   
#      of Technology, KanGAL Report 2001001                                

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
AMALGAMPar['N'] = 250    # Population size
AMALGAMPar['T'] = 500    # Number of generations
AMALGAMPar['d'] = 10     # Number of parameters
AMALGAMPar['m'] = 3      # Number of objective functions

# Parameter information for Latin hypercube sampling-based initialization
Par_info = {}
Par_info['initial'] = 'latin'                   # Latin hypercube sampling
Par_info['boundhandling'] = 'bound'             # Boundary handling
Par_info['min'] = np.zeros(AMALGAMPar['d'])     # Minimum parameter values
Par_info['max'] = np.ones(AMALGAMPar['d'])      # Maximum parameter values

# Define the function name
Func_name = 'AMALGAM_DTLZ2.AMALGAM_DTLZ2'       # The benchmark function

# Plugin for the second argument of Func_name (number of objectives)
plugin = AMALGAMPar['m']

# Set optimization options
options = {}
options['print'] = 'yes'   # Print output to screen (figures)
options['parallel'] = 'no'

# Compute the true Pareto front for the DTLZ2 problem
FX_true = Compute_FX_true(AMALGAMPar, Func_name, 'DTLZ2', 500, plugin)

if __name__ == '__main__':
    # Run the AMALGAM optimization and obtain non-dominated solution set
    X, FX, output, Z, YX = AMALGAM(AMALGAMPar, Func_name, Par_info, options, plugin, FX_true)
