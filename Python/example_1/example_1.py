# ####################################################################### #
#                                                                         #
#   EEEEEE  XX  XX   AAAA   MM   MM  PPPPPP  LL      EEEEEE        1111   #
#   EE       XXXX   AA  AA  MMM MMM  PP  PP  LL      EE           11 11   #
#   EEEEE     XX    AA  AA  MMMMMMM  PPPPPP  LL      EEEEE       11  11   #
#   EE       XXXX   AAAAAA  MM   MM  PP      LL      EE              11   #
#   EEEEEE  XX  XX  AA  AA  MM   MM  PP      LLLLLL  EEEEEE          11   #
#                                                                         #
# ####################################################################### #

# ZDT1 example from the following paper
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
AMALGAMPar = {'N': 100,    # Population size
			  'T': 100,    # Number of generations
			  'd': 30,     # Number of parameters
              'm': 2}      # Number of objective functions

# Parameter information for Latin Hypercube sampling
Par_info = {}
Par_info['initial'] = 'latin'           		# Latin hypercube sampling
Par_info['boundhandling'] = 'bound'    			# Boundary handling
Par_info['min'] = np.zeros(AMALGAMPar['d'])  	# Min values
Par_info['max'] = np.ones(AMALGAMPar['d'])  	# Max values

# Define function name
Func_name = 'AMALGAM_ZDT1.AMALGAM_ZDT1'  		# The benchmark function

# Plugin for the second argument of Func_name (number of parameters)
plugin = AMALGAMPar['d']

# Compute the true Pareto front for the ZDT1 problem
FX_true = Compute_FX_true(AMALGAMPar, Func_name, 'ZDT1', 500, plugin)

# Define computational settings
options = {'parallel': 'yes',			# Do not evaluate population in parallel
			'print': 'yes', 			# Print output to screen (figures)
            'IO': 'no',					# No input/output writing used by model
			'modout': 'yes'}

if __name__ == '__main__':
	# Run the AMALGAM optimization and obtain non-dominated solution set
	X, FX, output, Z, YX = AMALGAM(AMALGAMPar, Func_name, Par_info, options, plugin, FX_true)
