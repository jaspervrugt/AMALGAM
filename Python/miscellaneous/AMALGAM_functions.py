# ####################################################################### #
#                                                                         #
#     AAA    MMM    MMM    AAA    LLL      GGGGGGGGG    AAA    MMM    MMM #
#    AA AA   MMM    MMM   AA AA   LLL      GGGGGGGGG   AA AA   MMM    MMM #
#   AAA AAA  MMM    MMM  AAA AAA  LLL      GGG   GGG  AAA AAA  MMM    MMM #
#  AAA   AAA MMMM  MMMM AAA   AAA LLL      GGG   GGG AAA   AAA MMMM  MMMM #
#  AAA   AAA MMMMMMMMMM AAA   AAA LLL      GGGGGGGGG AAA   AAA MMMMMMMMMM #
#  AAAAAAAAA MMMMMMMMMM AAAAAAAAA LLL      GGGGGGGGG AAAAAAAAA MMMMMMMMMM #
#  AAAAAAAAA MMM    MMM AAAAAAAAA LLL            GGG AAAAAAAAA MMM    MMM #
#  AAA   AAA MMM    MMM AAA   AAA LLL            GGG AAA   AAA MMM    MMM #
#  AAA   AAA MMM    MMM AAA   AAA LLLLLLLL       GGG AAA   AAA MMM    MMM #
#  AAA   AAA MMM    MMM AAA   AAA LLLLLLLL GGGGGGGGG AAA   AAA MMM    MMM #
#                                                                         #
# ####################################################################### #
#                                                                         #
# AMALGAM: This general purpose MATLAB code is designed to find parameter #
# values that defines the Pareto trade-off surface corresponding to a     #
# vector of different objective functions. In principle, each Pareto      #
# solution is a different weighting of the objectives used. Therefore,    #
# one could use multiple trials with a single objective optimization      #
# algorithms using diferent values of the weights to find different       #
# Pareto solutions. However, various contributions to the optimization    #
# literature have demonstrated that this approach is rather inefficient.  #
# The AMALGAM code developed herein is designed to find an approximation  #
# of the Pareto solution set within a single optimization run. The        #
# AMALGAM method combines two new concepts, simultaneous multimethod      #
# search, and self-adaptive offspring creation, to ensure a fast,         #
# reliable, and computationally efficient solution to multiobjective      #
# optimization problems. This method is called a multi-algorithm,         #
# genetically adaptive multiobjective, or AMALGAM, method, to evoke the   #
# image of a procedure that blends the attributes of the best available   #
# individual optimization algorithms                                      #
#                                                                         #
# ####################################################################### #
#                                                                         #
#  ALGORITHM HAS BEEN DESCRIBED IN                                        #
#   Vrugt, J.A., Multi-criteria optimization using the AMALGAM software   #
#       package: Theory, concepts, and MATLAB implementation, UCI, 2015   #
#   Vrugt, J.A., B.A. Robinson, and J.M. Hyman (2009), Self-adaptive      #
#       multimethod search for global optimization in real-parameter      #
#       spaces, IEEE Transactions on Evolutionary Computation, 13(2),     #
#       pp. 243-259, https://doi.org/10.1109/TEVC.2008.924428             #
#   Vrugt, J.A., and B.A. Robinson (2007), Improved evolutionary          #
#       optimization from genetically adaptive multimethod search,        #
#       Proceedings of the National Academy of Sciences of the United     #
#       States of America, 104, pp. 708-711,                              #
#       https://doi.org/10.1073/pnas.061047110407                         #
#  MORE INFORMATION IN                                                    #
#   Vrugt, J.A., H.V. Gupta, L.A. Bastidas, W. Bouten, and S. Sorooshian  #
#       (2003), Effective and efficient algorithm for multi-objective     #
#       optimization of hydrologic models, Water Resources Research,      #
#       39(8), art. No. 1214, https://doi.org/10.1029/2002WR001746        #
#   Schoups, G.H., J.W. Hopmans, C.A. Young, J.A. Vrugt, and              #
#       W.W. Wallender, Multi-objective optimization of a regional        #
#       spatially-distributed subsurface water flow model, Journal of     #
#       Hydrology, pp. 20-48, 311(1-4),                                   #
#       https://doi.org/10.1016/j.jhydrol.2005.01.001, 2005.              #
#   Vrugt, J.A., P.H. Stauffer, T. Wöhling, B.A. Robinson, and            #
#       V.V. Vesselinov (2008), Inverse modeling of subsurface flow and   #
#       transport properties: A review with new developments, Vadose      #
#       Zone Journal, 7(2), 843 - 864,                                    #
#       https://doi.org/10.2136/vzj2007.0078                              #
#   Wöhling, T., J.A. Vrugt, and G.F. Barkle (2008), Comparison of three  #
#       multiobjective optimization algorithms for inverse modeling of    #
#       vadose zone hydraulic properties, Soil Science Society of America #
#       Journal, 72, 305 - 319, https://doi.org/10.2136/sssaj2007.0176    #
#   Wöhling, T., and J.A. Vrugt (2008), Combining multi-objective         #
#       optimization and Bayesian model averaging to calibrate forecast   #
#       ensembles of soil hydraulic models, Water Resources Research, 44, #
#       W12432, https://doi.org/10.1029/2008WR007154                      #
#                                                                         #
# ######################################################################  #
#                                                                         #
#  BUILT-IN CASE STUDIES                                                  #
#   Example 1   Multivariate normal benchmark study                       #
#   Example 1   ZDT1: test function                                       #
#   Example 2   ZDT2: test function                                       #
#   Example 3   ZDT3: test function                                       #
#   Example 4   ZDT4: test function                                       #
#   Example 5   ZDT6: test function                                       #
#   Example 6   ZDT6: test function, discrete parameter space             #
#   Example 7   DTLZ1: test function, 3 objectives                        #
#   Example 8   DTLZ2: test function, 3 objectives                        #
#   Example 9   DTLZ3: test function, 3 objectives                        #
#   Example 11  Real-world example using rainfall-discharge modeling      #
#   Example 12  Watershed modeling using driven & nondriven hydrograph    #
#   Example 13  Bayesian model averaging: RMSE, IS and CRPS               #
#   Example 14  Multi-criteria BMA training temperature ensemble          #
#   Example 15  Multi-criteria BMA training sea-level pressure ensemble   #
#                                                                         #
# ####################################################################### #
#                                                                         #
# COPYRIGHT (c) 2024  the author                                          #
#                                                                         #
#   This program is free software: you can modify it under the terms      #
#   of the GNU General Public License as published by the Free Software   #
#   Foundation, either version 3 of the License, or (at your option)      #
#   any later version                                                     #
#                                                                         #
#   This program is distributed in the hope that it will be useful, but   #
#   WITHOUT ANY WARRANTY; without even the implied warranty of            #
#   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU      #
#   General Public License for more details                               #
#                                                                         #
# ####################################################################### #
#                                                                         #
#  PYTHON CODE:                                                           #
#  © Written by Jasper A. Vrugt using GPT-4 OpenAI's language model       # 
#    University of California Irvine                                      #
#  Version 2.0    Feb 2025                                                #
#                                                                         #
# ####################################################################### #

import numpy as np
import pandas as pd
import seaborn as sns
import array, os, warnings, random, shutil, platform
import matplotlib.pyplot as plt
import matplotlib.animation as animation
from datetime import datetime
from itertools import combinations
from matplotlib.backends.backend_pdf import PdfPages
from screeninfo import get_monitors
from matplotlib.ticker import MaxNLocator
from matplotlib.lines import Line2D
import multiprocess as mp
import importlib

def AMALGAM_check(Func_name, AMALGAMPar, Par_info, options, Ftrue):
    """
    Check for setup errors.
    
    Parameters:
    Func_name (str): Name of the function or script.
    AMALGAMPar (dict): Parameters for the AMALGAM algorithm.
    Par_info (dict): Parameter information.
    options (dict): Options for the AMALGAM algorithm.
    Ftrue: Unused in the provided code but included for compatibility.
    
    Returns:
    tuple: Updated AMALGAMPar, Par_info, options.
    """
    # Derive current time and set deadline
    deadline = datetime.strptime('28-Feb-2025', '%d-%b-%Y')
    
    # Check if the trial version has expired
    if (deadline - datetime.now()).days < 0:
        raise ValueError('AMALGAM ERROR: Trial version of AMALGAM V2.0 has ended')
    
    # Open an output file with warnings
    with open('warning_file.txt', 'w+') as fid:
        fid.write('-------------- AMALGAM warning file --------------\n')

        # Check Func_name
        if not Func_name:
            raise ValueError("AMALGAM ERROR: Input argument Func_name should not be empty but a string enclosed between quotes with the name of the MATLAB function that calculates objective functions")
        
        if isinstance(Func_name, (int, float)):
            raise ValueError(f"AMALGAM ERROR: The variable Func_name is defined as a numerical value {Func_name}: This should be a string (between quotes) with the name of the MATLAB model script (.m file)")
        
        if not isinstance(Func_name, str):
            raise ValueError("AMALGAM ERROR: Input argument Func_name should be a string enclosed between quotes with the name of the MATLAB function that calculates objective functions")

        # Check AMALGAMPar
        if not AMALGAMPar:
            raise ValueError("AMALGAM ERROR: Input argument AMALGAMPar should not be empty but a dictionary with different fields")
        
        if not isinstance(AMALGAMPar, dict):
            raise ValueError("AMALGAM ERROR: Input argument AMALGAMPar should be a dictionary with fields")

        # Check Par_info
        if not Par_info:
            raise ValueError("AMALGAM ERROR: Input argument Par_info should not be empty but a dictionary with different fields")
        
        if not isinstance(Par_info, dict):
            raise ValueError("AMALGAM ERROR: Input argument Par_info should be a dictionary with fields")
        
        # Check options
        if not options:
            warnings.warn("AMALGAM WARNING: Input argument options is empty -> default settings will be assumed for its fields (see Table 3 of the manual)")
            fid.write("AMALGAM WARNING: Input argument options is empty -> default settings will be assumed for its fields (see Table 3 of manual)\n")
        
        if not isinstance(options, dict):
            raise ValueError("AMALGAM ERROR: Input argument options should be a dictionary with fields")
        
        # Convert all strings to lowercase except parameter names
        for structure in ['AMALGAMPar', 'Par_info', 'options']:
            vrbl = eval(structure)
            for key, value in vrbl.items():
                if isinstance(value, str):
                    vrbl[key] = value.lower()
    
        # Validate 'N' (population size)
        if 'N' not in AMALGAMPar:
            raise ValueError("AMALGAM ERROR: Field 'N' of structure AMALGAMPar undefined")
        elif not AMALGAMPar['N']:
            raise ValueError("AMALGAM ERROR: Field 'N' of structure AMALGAMPar left empty - please list population size")
        elif not isinstance(AMALGAMPar['N'], int):
            raise ValueError("AMALGAM WARNING: Field 'N' of structure AMALGAMPar should be an integer - please list population size")

        AMALGAMPar['N'] = int(AMALGAMPar['N'])

        # Validate 'd' (dimensionality)
        if 'd' not in AMALGAMPar:
            raise ValueError("AMALGAM ERROR: Field 'd' of structure AMALGAMPar undefined")
        elif not AMALGAMPar['d']:
            raise ValueError("AMALGAM ERROR: Field 'd' of structure AMALGAMPar left empty - please list problem dimensionality")
        elif not isinstance(AMALGAMPar['d'], int):
            raise ValueError("AMALGAM ERROR: Field 'd' of structure AMALGAMPar should be an integer - please list problem dimensionality")

        AMALGAMPar['d'] = int(AMALGAMPar['d'])

        # Validate 'T' (maximum number of generations)
        if 'T' not in AMALGAMPar:
            raise ValueError("AMALGAM ERROR: Field 'T' of structure AMALGAMPar undefined")
        elif not AMALGAMPar['T']:
            raise ValueError("AMALGAM ERROR: Field 'T' of structure AMALGAMPar left empty - please list the maximum number of generations")
        elif not isinstance(AMALGAMPar['T'], int):
#            raise ValueError("AMALGAM ERROR: Field 'T' of structure AMALGAMPar should be an integer - please list maximum number of generations")
#            print("AMALGAM WARNING: Field 'T' of structure AMALGAMPar should be an integer - please list maximum number of generations")
#            raise ValueError("AMALGAM ERROR: Field 'd' of structure AMALGAMPar should be an integer - please list problem dimensionality")
            AMALGAMPar['T'] = int(AMALGAMPar['T'])  
        # Validate 'm' (number of objective functions)
        if 'm' not in AMALGAMPar:
            raise ValueError("AMALGAM ERROR: Field 'm' of structure AMALGAMPar undefined")
        elif not AMALGAMPar['m']:
            raise ValueError("AMALGAM ERROR: Field 'm' of structure AMALGAMPar left empty - please list the number of objective functions")
        elif not isinstance(AMALGAMPar['m'], int):
#            raise ValueError("AMALGAM ERROR: Field 'm' of structure AMALGAMPar should be an integer - please list the number of objective functions")
#            print("AMALGAM WARNING: Field 'm' of structure AMALGAMPar should be an integer - please list the number of objective functions")
#            raise ValueError("AMALGAM ERROR: Field 'd' of structure AMALGAMPar should be an integer - please list problem dimensionality")
            AMALGAMPar['m'] = int(AMALGAMPar['m'])
        # Check population size
        if AMALGAMPar['N'] < 30:
            raise ValueError("AMALGAM ERROR: Recommend to use a larger population size -> Set AMALGAMPar.N to be at least 50")
        # Check dimensionality
        if AMALGAMPar['d'] <= 0:
            raise ValueError("AMALGAM ERROR: Number of parameters should be integer and larger than zero -> Set AMALGAMPar['d'] >= 1")
        # Check generations
        if AMALGAMPar['T'] < 2:
            raise ValueError("AMALGAM ERROR: Number of generations smaller than one -> Set at least AMALGAMPar.T = 2")
        # Validate 'rec_methods' (recombination methods)
        if 'rec_methods' in AMALGAMPar:
            if not AMALGAMPar['rec_methods']:
                print("AMALGAM WARNING: Field 'rec_methods' of structure AMALGAMPar left empty -> default set of recombination methods is used")
            elif not isinstance(AMALGAMPar['rec_methods'], list):
                raise ValueError("AMALGAM ERROR: Field 'rec_methods' of structure AMALGAMPar should be a list with the names (acronyms) of the recombination methods used")
            else:
                valid_rec_methods = ['ga', 'ps', 'am', 'de']
                for rec_method in AMALGAMPar['rec_methods']:
                    if rec_method not in valid_rec_methods:
                        raise ValueError(f"AMALGAM ERROR: Unknown recombination method {rec_method} --> Select from {valid_rec_methods}")
        # Validate 'beta_1' (DE scaling factor)
        if 'beta_1' in AMALGAMPar:
            if not AMALGAMPar['beta_1']:
                print("AMALGAM WARNING: Field 'beta_1' of structure AMALGAMPar left empty - default value assumed")
            elif not isinstance(AMALGAMPar['beta_1'], (int, float)):
                raise ValueError("AMALGAM ERROR: Field 'beta_1' of structure AMALGAMPar should be a numerical value (default: AMALGAMPar.beta_1 drawn from UNIFORM[0.6,1])")
            if AMALGAMPar['beta_1'] <= 0:
                raise ValueError("AMALGAM ERROR: DE scaling factor beta_1 of differential evolution should be larger than zero -> Use at least AMALGAMPar.beta_1 = 0.4")
            elif AMALGAMPar['beta_1'] > 2:
                print("AMALGAM WARNING: DE scaling factor beta_1 of differential evolution somewhat large --> (default: AMALGAMPar.beta_1 drawn from UNIFORM[0.6,1])")   
        # Validate 'beta_2' (DE scaling factor)
        if 'beta_2' in AMALGAMPar:
            if not AMALGAMPar['beta_2']:
                print("AMALGAM WARNING: Field 'beta_2' of structure AMALGAMPar left empty - default value assumed")
            elif not isinstance(AMALGAMPar['beta_2'], (int, float)):
                raise ValueError("AMALGAM ERROR: Field 'beta_2' of structure AMALGAMPar should be a numerical value (default: AMALGAMPar.beta_2 drawn from UNIFORM[0.2,0.6])")
            if AMALGAMPar['beta_2'] <= 0:
                raise ValueError("AMALGAM ERROR: DE scaling factor beta_2 of differential evolution should be larger than zero -> Use at least AMALGAMPar.beta_2 = 0.2")
            elif AMALGAMPar['beta_2'] > 2:
                print("AMALGAM WARNING: DE scaling factor beta_2 of differential evolution set rather large --> (default: AMALGAMPar.beta_2 drawn from UNIFORM[0.2,0.6])")
        # Validate 'c1'
        if 'c1' in AMALGAMPar:
            if AMALGAMPar['c1'] is None:
                evalstr = "AMALGAM WARNING: Field 'c1' of structure AMALGAMPar left empty - default value assumed (see Table 1 of manual)\n"
                print(evalstr)
                if fid:
                    fid.write(evalstr)
            elif not isinstance(AMALGAMPar['c1'], (int, float)):
                raise ValueError("AMALGAM ERROR: Field 'c1' of structure AMALGAMPar should be a numerical value (default: AMALGAMPar.c1 = 1.5)")
        
            if AMALGAMPar['c1'] <= 0:
                raise ValueError("AMALGAM ERROR: PSO social factor should be larger than zero -> Use at least AMALGAMPar.c1 = 1 (default: AMALGAMPar.c1 = 1.5)")
            elif AMALGAMPar['c1'] > 2:
                evalstr = "AMALGAM WARNING: PSO social factor set rather large -> (default: AMALGAMPar.c1 = 1.5)\n"
                print(evalstr)
            if fid:
                fid.write(evalstr)
        # Validate 'c2'
        if 'c2' in AMALGAMPar:
            if AMALGAMPar['c2'] is None:
                evalstr = "AMALGAM WARNING: Field 'c2' of structure AMALGAMPar left empty - default value assumed (see Table 1 of manual)\n"
                print(evalstr)
                if fid:
                    fid.write(evalstr)
            elif not isinstance(AMALGAMPar['c2'], (int, float)):
                raise ValueError("AMALGAM ERROR: Field 'c2' of structure AMALGAMPar should be a numerical value (default: AMALGAMPar.c2 = 1.5)")

            if AMALGAMPar['c2'] <= 0:
                raise ValueError("AMALGAM ERROR: PSO cognitive factor should be larger than zero -> Use at least AMALGAMPar.c2 = 1 (default: AMALGAMPar.c2 = 1.5)")
            elif AMALGAMPar['c2'] > 2:
                evalstr = "AMALGAM WARNING: PSO cognitive factor set rather large -> (default: AMALGAMPar.c2 = 1.5)\n"
                print(evalstr)
                if fid:
                    fid.write(evalstr)
        # Validate 'varphi'
        if 'varphi' in AMALGAMPar:
            if AMALGAMPar['varphi'] is None:
                evalstr = "AMALGAM WARNING: Field 'varphi' of structure AMALGAMPar left empty - default value assumed (see Table 1 of manual)\n"
                print(evalstr)
                if fid:
                    fid.write(evalstr)
            elif not isinstance(AMALGAMPar['varphi'], (int, float)):
                raise ValueError("AMALGAM ERROR: Field 'varphi' of structure AMALGAMPar should be a numerical value between 0 and 1 (default: AMALGAMPar.varphi drawn from UNIFORM[0.5,1.0])")

            if AMALGAMPar['varphi'] <= 0:
                raise ValueError("AMALGAM ERROR: PSO inertia factor should be larger than zero -> Use at least AMALGAMPar.varphi = 0.25 (default: AMALGAMPar.varphi drawn from UNIFORM[0.5,1.0])")
            elif AMALGAMPar['varphi'] > 2:
                evalstr = "AMALGAM WARNING: PSO inertia factor set rather large -> (default: AMALGAMPar.varphi drawn from UNIFORM[0.5,1.0])\n"
                print(evalstr)
                if fid:
                    fid.write(evalstr)
        # Validate 'p_CR'
        if 'p_CR' in AMALGAMPar:
            if AMALGAMPar['p_CR'] is None:
                evalstr = "AMALGAM WARNING: Field 'p_CR' of structure AMALGAMPar left empty - default value assumed (see Table 1 of manual)\n"
                print(evalstr)
                if fid:
                    fid.write(evalstr)
            elif not isinstance(AMALGAMPar['p_CR'], (int, float)):
                raise ValueError("AMALGAM ERROR: Field 'p_CR' of structure AMALGAMPar should be a numerical value between 0 and 1 (default: AMALGAMPar.p_CR = 0.9)")

            if AMALGAMPar['p_CR'] < 0 or AMALGAMPar['p_CR'] > 1:
                raise ValueError("AMALGAM ERROR: NSGA-II crossover rate should be between 0 and 1 -> AMALGAMPar.p_CR in [0,1] (default: AMALGAMPar.p_CR = 0.9)")
        # Validate 'p_M'
        if 'p_M' in AMALGAMPar:
            if AMALGAMPar['p_M'] is None:
                evalstr = "AMALGAM WARNING: Field 'p_M' of structure AMALGAMPar left empty - default value assumed (see Table 1 of manual)\n"
                print(evalstr)
                if fid:
                    fid.write(evalstr)
            elif not isinstance(AMALGAMPar['p_M'], (int, float)):
                raise ValueError("AMALGAM ERROR: Field 'p_M' of structure AMALGAMPar should be a numerical value (default: AMALGAMPar.p_M = 1/AMALGAMPar['d'])")

            if AMALGAMPar['p_M'] < 0 or AMALGAMPar['p_M'] > 1:
                raise ValueError("AMALGAM ERROR: NSGA-II mutation rate should be between 0 and 1 -> AMALGAMPar.p_M in [0,1] (default: AMALGAMPar.p_M = 1/AMALGAMPar['d'])")
        # Validate 'eta_C'
        if 'eta_C' in AMALGAMPar:
            if AMALGAMPar['eta_C'] is None:
                evalstr = "AMALGAM WARNING: Field 'eta_C' of structure AMALGAMPar left empty - default value assumed (see Table 1 of manual)\n"
                print(evalstr)
                if fid:
                    fid.write(evalstr)
            elif not isinstance(AMALGAMPar['eta_C'], (int, float)):
                raise ValueError("AMALGAM ERROR: Field 'eta_C' of structure AMALGAMPar should be a numerical value (default: AMALGAMPar.eta_C = 10)")

            if AMALGAMPar['eta_C'] <= 0:
                raise ValueError("AMALGAM ERROR: NSGA-II mutation index eta_C should be larger than zero -> Use at least AMALGAMPar.eta_C = 1 (default: AMALGAMPar.eta_C = 10)")
            elif AMALGAMPar['eta_C'] > 200:
                evalstr = "AMALGAM WARNING: NSGA-II mutation index eta_C set rather large -> (default: AMALGAMPar.eta_C = 10)\n"
                print(evalstr)
                if fid:
                    fid.write(evalstr)
        # Validate 'eta_M'
        if 'eta_M' in AMALGAMPar:
            if AMALGAMPar['eta_M'] is None:
                evalstr = "AMALGAM WARNING: Field 'eta_M' of structure AMALGAMPar left empty - default value assumed (see Table 1 of manual)\n"
                print(evalstr)
                if fid:
                    fid.write(evalstr)
            elif not isinstance(AMALGAMPar['eta_M'], (int, float)):
                raise ValueError("AMALGAM ERROR: Field 'eta_M' of structure AMALGAMPar should be a numerical value (default: AMALGAMPar.eta_M = 50)")

            if AMALGAMPar['eta_M'] <= 0:
                raise ValueError("AMALGAM ERROR: NSGA-II mutation index eta_M should be larger than zero -> Use at least AMALGAMPar.eta_M = 1 (default: AMALGAMPar.eta_M = 50)")
            elif AMALGAMPar['eta_M'] > 200:
                evalstr = "AMALGAM WARNING: PSO mutation index eta_M set rather large -> (default: AMALGAMPar.eta_M = 50)\n"
                print(evalstr)
                if fid:
                    fid.write(evalstr)
        # Validate 'gamma'
        if 'gamma' in AMALGAMPar:
            if AMALGAMPar['gamma'] is None:
                evalstr = "AMALGAM WARNING: Field 'gamma' of structure AMALGAMPar left empty - default value assumed (see Table 1 of manual)\n"
                print(evalstr)
                if fid:
                    fid.write(evalstr)
            elif not isinstance(AMALGAMPar['gamma'], (int, float)):
                raise ValueError("AMALGAM ERROR: Field 'gamma' of structure AMALGAMPar should be a numerical value (default: AMALGAMPar.gamma = (2.38/sqrt(AMALGAMPar['d']))^2)")

            if AMALGAMPar['gamma'] <= 0:
                raise ValueError("AMALGAM ERROR: AMS jump rate should be larger than zero -> Use at least AMALGAMPar.gamma = 0.1 (default: AMALGAMPar.gamma = (2.38/sqrt(AMALGAMPar['d']))^2)")
            elif AMALGAMPar['gamma'] > 5:
                evalstr = "AMALGAM WARNING: AMS jump rate set rather large -> (default: AMALGAMPar.gamma = (2.38/sqrt(AMALGAMPar['d']))^2)\n"
                print(evalstr)
                if fid:
                    fid.write(evalstr)
        # Validate 'K'
        if 'K' in AMALGAMPar:
            if AMALGAMPar['K'] is None:
                evalstr = "AMALGAM WARNING: Field 'K' of structure AMALGAMPar left empty - default value assumed (see Table 1 of manual)\n"
                print(evalstr)
                if fid:
                    fid.write(evalstr)
            elif not isinstance(AMALGAMPar['K'], int):
                raise ValueError("AMALGAM ERROR: Field 'K' of structure AMALGAMPar should be an integer (default: K = 1)")

            if AMALGAMPar['K'] < 1:
                raise ValueError("AMALGAM ERROR: Thinning parameter should be integer and larger than zero -> Set AMALGAMPar.K >= 1 (default: 1)")
            elif AMALGAMPar['K'] > 20:
                evalstr = "AMALGAM WARNING: Field 'K' of structure AMALGAMPar set rather large -> recommend to use AMALGAMPar.K in [1,20]\n"
                print(evalstr)
                if fid:
                    fid.write(evalstr)
        # Validate 'p0'
        if 'p0' in AMALGAMPar:
            if AMALGAMPar['p0'] is None:
                evalstr = "AMALGAM WARNING: Field 'p0' of structure AMALGAMPar left empty - default value assumed (see Table 1 of manual)\n"
                print(evalstr)
                if fid:
                    fid.write(evalstr)
            elif not isinstance(AMALGAMPar['p0'], (int, float)):
                raise ValueError("AMALGAM ERROR: Field 'p0' of structure AMALGAMPar should contain a numerical value")

            if not np.isscalar(AMALGAMPar['p0']):
                raise ValueError("AMALGAM ERROR: Too many minimum values of selection probabilities. AMALGAMPar.p0 should be a scalar")

            if AMALGAMPar['p0'] < 0:
                raise ValueError("AMALGAM ERROR: Minimum selection probability of recombination methods should be larger than zero")

            if 'rec_methods' in AMALGAMPar:
                n_scalar = 1 / len(AMALGAMPar['rec_methods'])
            else:
                n_scalar = 1 / 4

            if AMALGAMPar['p0'] > n_scalar:
                raise ValueError(f"AMALGAM ERROR: Minimum selection probability of recombination methods set too large --> AMALGAMPar.p0 in [0, {round(100 * n_scalar) / 100}]")

            if AMALGAMPar['p0'] == n_scalar:
                evalstr = f"AMALGAM WARNING: Minimum selection probability set so large that no adaptation can take place --> AMALGAMPar.p0 in [0, {round(100 * n_scalar) / 200}] to get some adaptation"
                print(evalstr)
                if fid:
                    fid.write(evalstr)
        # Validate if 'initial' is defined
        if 'initial' not in Par_info:
            raise ValueError("AMALGAM ERROR: Initial sampling distribution not defined -> Define Par_info.initial = 'latin' or 'uniform' or 'normal' or 'prior'!!")
        elif not Par_info['initial']:
            raise ValueError("AMALGAM ERROR: Field 'initial' of structure Par_info left empty -> Define Par_info.initial = 'latin' or 'uniform' or 'normal' or 'prior'!!")
        elif not isinstance(Par_info['initial'], str):
            raise ValueError("AMALGAM ERROR: Field 'initial' of structure Par_info should contain a string enclosed between quotes -> Par_info.initial = 'latin' or 'uniform' or 'normal' or 'prior'!!")  
        # Validate if 'initial' is one of the valid options
        if Par_info['initial'] not in ['latin', 'uniform', 'normal', 'prior']:
            raise ValueError("AMALGAM ERROR: Initial sampling distribution unknown -> Set Par_info.initial = 'latin' or 'uniform' or 'normal' or 'prior'!!")
        # Latin hypercube sampling requires 'min' and 'max' fields
        if Par_info['initial'] == 'latin':
            if 'min' not in Par_info:
                raise ValueError("AMALGAM ERROR: Latin hypercube sampling selected but minimum parameter values not defined -> Set Par_info.min!!")
            if 'max' not in Par_info:
                raise ValueError("AMALGAM ERROR: Latin hypercube sampling selected but maximum parameter values not defined -> Set Par_info.max!!")
        # Uniform sampling requires 'min' and 'max' fields
        if Par_info['initial'] == 'uniform':
            if 'min' not in Par_info:
                raise ValueError("AMALGAM ERROR: Uniform initial sampling selected but minimum parameter values not defined -> Set Par_info.min!!")
            if 'max' not in Par_info:
                raise ValueError("AMALGAM ERROR: Uniform initial sampling selected but maximum parameter values not defined -> Set Par_info.max!!")
        # Normal distribution requires 'mu' and 'cov' fields
        if Par_info['initial'] == 'normal':
            if 'mu' not in Par_info:
                raise ValueError("AMALGAM ERROR: Normal distribution selected to sample from but unknown mean -> Define Par_info.mu!!")
            if 'cov' not in Par_info:
                raise ValueError("AMALGAM ERROR: Normal distribution selected to sample from but unknown covariance -> Define Par_info.cov!!")   
        # Validate 'mu' and 'cov' for normal distribution
        if Par_info['initial'] == 'normal':
            if len(Par_info['mu']) != AMALGAMPar['d']:
                raise ValueError(f"AMALGAM ERROR: Mean of normal distribution (Par_info.mu) should be a row vector with AMALGAMPar['d'] = {AMALGAMPar['d']} values")
            if len(Par_info['cov']) != AMALGAMPar['d'] or len(Par_info['cov'][0]) != AMALGAMPar['d']:
                raise ValueError(f"AMALGAM ERROR: Covariance of normal distribution ('Par_info.cov') should be a square matrix of size AMALGAMPar['d'] x AMALGAMPar['d'] = {AMALGAMPar['d']} x {AMALGAMPar['d']} values")
        # Prior distribution requires 'prior' field
        if Par_info['initial'] == 'prior':
            if 'prior' not in Par_info:
                raise ValueError("AMALGAM ERROR: Prior distribution selected but unknown field 'prior' of structure Par_info -> Define Par_info.prior!!")
        # Check for 'boundhandling'
        if 'boundhandling' not in Par_info:
            evalstr = "AMALGAM CHECK: Field 'boundhandling' of structure Par_info not specified by user -> thus AMALGAM assumes an unbounded parameter space\n"
            print(evalstr)
            fid.write(evalstr)
            Par_info['boundhandling'] = 'none'
        elif not Par_info['boundhandling']:
            evalstr = "AMALGAM WARNING: Field 'boundhandling' of structure Par_info is empty -> AMALGAM assumes unbounded parameter space\n"
            print(evalstr)
            fid.write(evalstr)
            Par_info['boundhandling'] = 'none'
        elif not isinstance(Par_info['boundhandling'], str):
            raise ValueError("AMALGAM ERROR: Field 'boundhandling' of structure Par_info should contain a string enclosed between quotes -> Par_info.boundhandling = 'fold' or 'bound' or 'reflect' or 'none'!!")

        if Par_info['boundhandling'] not in ['fold', 'bound', 'reflect', 'none']:
            raise ValueError("AMALGAM ERROR: Unknown boundary handling method -> Define Par_info.boundhandling = 'fold' or 'bound' or 'reflect' or 'none'!!")
        # If boundary handling is specified, 'min' and 'max' must be defined
        if Par_info['boundhandling'] in ['fold', 'bound', 'reflect']:
            if 'min' not in Par_info:
                raise ValueError(f"AMALGAM ERROR: Par_info.boundhandling is equal to {Par_info['boundhandling']} but minimum parameter values not defined in Par_info.min")
            if 'max' not in Par_info:
                raise ValueError(f"AMALGAM ERROR: Par_info.boundhandling is equal to {Par_info['boundhandling']} but maximum parameter values not defined in Par_info.max")                         
            # added
            if isinstance(Par_info['min'], array.array):
                Par_info['min'] = np.array(Par_info['min']).reshape(1, -1)
            elif len(Par_info['min']) == 0:
                raise ValueError(f"AMALGAM ERROR: Par_info.boundhandling is equal to {Par_info['boundhandling']} but Par_info.min is empty -> please define minimum values of each parameter")
            elif not isinstance(Par_info['min'][0], (int, float, np.int64, np.float64)):
                raise ValueError(f"AMALGAM ERROR: Par_info.boundhandling is equal to {Par_info['boundhandling']} but Par_info.min does not list numerical lower limits parameters")
            # added
            if isinstance(Par_info['max'], array.array):
                Par_info['max'] = np.array(Par_info['max']).reshape(1,-1)
            elif len(Par_info['max']) == 0:
                raise ValueError(f"AMALGAM ERROR: Par_info.boundhandling is equal to {Par_info['boundhandling']} but Par_info.max is empty -> please define maximum values of each parameter")
            elif not isinstance(Par_info['max'][0], (int, float, np.int64, np.float64)):
                raise ValueError(f"AMALGAM ERROR: Par_info.boundhandling is equal to {Par_info['boundhandling']} but Par_info.max does not list numerical upper limits parameters")
            
        # Ensure that 'min' and 'max' are consistent with AMALGAMPar['d']
        Par_info['min'] = np.array(Par_info['min']).reshape(1,-1)
        Par_info['max'] = np.array(Par_info['max']).reshape(1,-1)
        if 'min' in Par_info and Par_info['min'].shape[1] != AMALGAMPar['d']:
            raise ValueError(f"AMALGAM ERROR: Number of elements of field 'min' of structure Par_info should be equal to {AMALGAMPar['d']}!!")
        if 'max' in Par_info and Par_info['max'].shape[1] != AMALGAMPar['d']:
            raise ValueError(f"AMALGAM ERROR: Number of elements of field 'max' of structure Par_info should be equal to {AMALGAMPar['d']}!!")
        # Validate 'steps' field
        if 'steps' in Par_info:
            # if not Par_info['steps']:
            if Par_info['steps'] is None:
                raise ValueError("AMALGAM ERROR: Field 'steps' of structure Par_info left empty -> Define Par_info.steps to have AMALGAMPar['d'] elements!!")
            if any([round(x) != x for x in Par_info['steps']]):
                raise ValueError("AMALGAM ERROR: Field 'steps' of structure Par_info should contain integers only!!")
            if len(Par_info['steps']) == 1 and AMALGAMPar['d'] > 1:
                evalstr = "AMALGAM WARNING: Field 'steps' of structure Par_info lists a single value and will be copied AMALGAMPar['d'] times\n"
                print(evalstr)
                fid.write(evalstr)
                Par_info['steps'] = [Par_info['steps'][0]] * AMALGAMPar['d']
            elif len(Par_info['steps']) != AMALGAMPar['d']:
                raise ValueError(f"AMALGAM ERROR: Number of elements of field 'steps' of structure Par_info should be equal to {AMALGAMPar['d']}!!")
        # Remove 'boundhandling' if not used
        if Par_info.get('boundhandling') == 'none':
            del Par_info['boundhandling']
        # Validate 'density' field in options
        if 'density' in options:
            if options['density'] == '':
                raise ValueError("AMALGAM ERROR: Field 'density' of structure options is empty -> Set options.density = 'crowding' or 'strength'")
            elif not isinstance(options['density'], str):
                raise ValueError("AMALGAM ERROR: Field 'density' of structure options should contain string enclosed between quotes -> Set options.density = 'crowding' or 'strength'")
            elif options['density'] not in ['crowding', 'strength']:
                raise ValueError("AMALGAM ERROR: Unknown density method -> Set options.density = 'crowding' or 'strength' (default 'crowding')")
        else:
            evalstr = ("AMALGAM WARNING: Field 'density' of structure options not defined -> resort to default setting, options.density = 'crowding'\n")
            # Print warning to screen and file
            print(evalstr)
            fid.write(evalstr)
        # Validate 'print' field in options
        if 'print' not in options:
            evalstr = "AMALGAM WARNING: Field 'print' of structure options not defined as 'yes' or 'no' -> resort to default setting of options.print = 'yes'\n"
            # Print warning to screen and file
            print(evalstr)
            fid.write(evalstr)
        else:
            if options['print'] == '':
                raise ValueError("AMALGAM ERROR: Field 'print' of structure options should not be empty but should contain a string ('yes' or 'no')")
            elif not isinstance(options['print'], str):
                raise ValueError("AMALGAM ERROR: Field 'print' of structure options should be a string ('yes' or 'no')")
            elif options['print'] not in ['yes', 'no']:
                raise ValueError("AMALGAM ERROR: Field 'print' of structure options should equal 'yes' or 'no'")
        # Validate 'screen' field in options
        if 'screen' not in options or options['screen'] == '' or not isinstance(options['screen'], str) or options['screen'] not in ['yes', 'no']:
            options['screen'] = 'no'
        # Validate each field of structure options
        for field in options:
            F = options[field]
            if field not in ['ranking', 'density']:
                if F not in ['yes', 'no']:
                    raise ValueError(f"AMALGAM ERROR: Field '{field}' of structure options should be set equal to 'yes' or 'no'")

        # Validate 'Ftrue' for Pareto solutions
        if Ftrue is None:
            evalstr = "AMALGAM WARNING: User did not specify Pareto solutions (sixth input argument is empty) -> cannot compute convergence diagnostics\n"
            print(evalstr)
            fid.write(evalstr)
        elif not isinstance(Ftrue, (list, np.ndarray)):
            raise ValueError("AMALGAM ERROR: Input argument Fpar should be left empty or contain numerical values (for benchmark experiments on KNOWN Pareto solution sets)")
        if Ftrue is not None and len(Ftrue) > 0: 
            if Ftrue.shape[1] != AMALGAMPar['m']:
                raise ValueError(f"AMALGAM ERROR: Number of columns of Fpar should match number of objectives, that is, AMALGAMPar['m'] = {AMALGAMPar['m']} columns in present case")
            elif Ftrue.shape[0] < 5:
                raise ValueError("AMALGAM ERROR: Number of pareto solutions stored in Fpar is rather small -> increase to at least 100 to get reasonable convergence statistics")

    # Close the warning file
    fid.close()

    return AMALGAMPar, Par_info, options


def AMALGAM_setup(AMALGAMPar, Par_info, options):
    # ####################################################################### #
    #                                                                         #
    # Initializes the main variables used in AMALGAM                          #
    #                                                                         #
    #  SYNOPSIS                                                               #
    #   [AMALGAMPar,Par_info,options,T_start] = AMALGAM_setup( ...            #
    #       AMALGAMPar,Par_info,options)                                      #
    #  where                                                                  #
    #                                                                         #
    # ####################################################################### #
    
    # Random seed (equivalent to 'rng(1+round(100*rand),'twister'); in MATLAB)
    random.seed(1 + round(100 * random.random()))

    # Define T_start
    T_start = 1

    # Field names of AMALGAMPar and their default values
    f_names = ['beta_1', 'beta_2', 'c_1', 'c_2', 'varphi', 'p_CR', 'p_M', 'eta_C', 'eta_M', 'gamma', 'K', 'p0']
    value = [
        #lambda N: np.random.uniform(0.6, 1.0, N).reshape(-1,1),
        0.8, # function is now local in DE for restart run: U[beta_1 - 0.2, beta_1 + 0.2]
        #lambda N: np.random.uniform(0.2, 0.6, N).reshape(-1,1),
        0.4, # function is now local in DE for restart run: U[beta_2 - 0.2, beta_2 + 0.2]
        1.5, 1.5,
        #lambda N: np.random.uniform(0.5, 1.0, N).reshape(-1,1),
        0.75, # function is now local in PSO for restart run: U[varphi - 1/4, varphi + 1/4] 
        0.9, 
        1 / AMALGAMPar['d'], 
        20, 
        20, 
        (2.38 / np.sqrt(AMALGAMPar['d'])) ** 2, 
        1, 
        0.05]
    
    # Check each field and set default value if missing
    for i, fname in enumerate(f_names):
        if fname not in AMALGAMPar:
            AMALGAMPar[fname] = value[i]

     # If recombination methods not specified --> use default AMALGAM
    if 'rec_methods' not in AMALGAMPar:
        AMALGAMPar['rec_methods'] = ['ga', 'ps', 'am', 'de']

    # Field names of options and their default values
    f_names = ['parallel', 'IO', 'save', 'restart', 'modout', 'density', 'print']
    value = ['no', 'no', 'no', 'no', 'no', 'crowding', 'yes']
    
    # Set undefined fields to default values
    for i, fname in enumerate(f_names):
        if fname not in options:
            options[fname] = value[i]

    # Replicate Par_info.min and Par_info.max (assuming they are numpy arrays)
    # Par_info['min'] = np.tile(Par_info['min'], (AMALGAMPar['N'], 1)) 
    # Par_info['max'] = np.tile(Par_info['max'], (AMALGAMPar['N'], 1))

    # Compute step size and make N copies
    if 'steps' in Par_info:
        # Par_info['steps'] = np.tile(Par_info['steps'], (AMALGAMPar['N'], 1))
        Par_info['step_size'] = (Par_info['max'] - Par_info['min']) / Par_info['steps']

    # Define prior handle (random samples & evaluation of pdf)
    if 'prior' in Par_info:
        if isinstance(Par_info['prior'], list):
            Par_info['u'] = 'yes'  # Univariate case
            Par_info['prior_rnd'] = []
            for ii in range(AMALGAMPar['d']):
                # Handle draw initial state
                prior_rnd_func = eval(f"lambda x: {Par_info['prior'][ii].replace('pdf', 'rnd')}")
                Par_info['prior_rnd'].append(prior_rnd_func)
        else:
            Par_info['u'] = 'no'  # Multivariate case
            pr_name = str(Par_info['prior'])
            pr_var, idcp = extract_names(pr_name)  # Assumes extract_names is defined elsewhere
            pr_name = pr_name[idcp[0] + 1:]  # Prior handle without @()
            n_var = len(pr_var)  # Number of variables without x
            for z in range(n_var):
                pr_name = pr_name.replace(pr_var[z], f"Par_info.{pr_var[z]}")
            Par_info['prior_rnd'] = eval(f"lambda x: {pr_name.replace('pdf(x', 'rnd(')}")

    # Check for BMA (assuming global BMA is defined elsewhere)
    if 'BMA' in options:
        if options['BMA'] == 'yes':
            Par_info['unit_simplex'] = 'yes'
    else:
        pass  

    # Order fields of AMALGAMPar and print summary
    AMALGAMPar = {k: AMALGAMPar[k] for k in ['d', 'N', 'T', 'm', 'rec_methods', 'K', 'p0', 'beta_1', 'beta_2', 'c_1', 'c_2', 'varphi', 'p_CR', 'p_M', 'eta_C', 'eta_M', 'gamma']}

    # Print summary of the settings if flag is 0
    flag = 0
    if flag == 0:
        print("---------- Summary of the main settings used: AMALGAMPar -----------------")
        for key, value in AMALGAMPar.items():
            if callable(value):
                print(f"{key:<11}: {value}")        
            elif isinstance(value, np.ndarray):                                   
                flattened_values = value.flatten()  
                formatted_values = ', '.join([f"{v:.4f}" for v in flattened_values]) 
                print(f"{key:<11}: {formatted_values}")
            elif isinstance(value, list):
                pr_str = '{' + ', '.join(map(str, value)) + '}'
                print(f"{key:<11}: {pr_str}")
            else:
                print(f"{key:<11}: {value:<6}")
        print("--------------------------------------------------------------------------")

        print("----------- Summary of the main settings used: Par_info ------------------")
        for key, value in Par_info.items():
            if callable(value):
                print(f"{key:<13}: {value}")
            elif isinstance(value, np.ndarray):  # Handle numpy arrays
                flattened_values = value.flatten()
                formatted_values = ', '.join([f"{v:.2f}" for v in flattened_values if isinstance(v, (int, float))])  # Only format numeric values
                print(f"{key:<13}: {formatted_values}")
            elif isinstance(value, list):  # Handle lists
                pr_str = '{' + ', '.join(map(str, value)) + '}'
                print(f"{key:<13}: {pr_str}")
            elif isinstance(value, (int, float)):  # Handle numeric types (int, float)
                print(f"{key:<13}: {value:.2f}")  # Apply numeric formatting only for numbers
            else:
                print(f"{key:<13}: {value}")  # Just print the value for other types                
        print("--------------------------------------------------------------------------")

        print("----------- Summary of the main settings used: options -------------------")
        for key, value in options.items():
            if callable(value):
                print(f"{key:<11}: {value}")        
            elif isinstance(value, np.ndarray):  
                flattened_values = value.flatten() 
                formatted_values = ', '.join([f"{v:.4f}" for v in flattened_values])  
                print(f"{key:<11}: {formatted_values}")
            elif isinstance(value, list):
                pr_str = '{' + ', '.join(map(str, value)) + '}'
                print(f"{key:<11}: {pr_str}")
            else:
                print(f"{key:<11}: {value:<6}")
        print("--------------------------------------------------------------------------")
    
    # Print header information
    print("\n")
    print('  -----------------------------------------------------------------------------------                ')
    print('      AAA     MMM        MMM     AAA     LLL     GGGGGGGGG     AAA     MMM        MMM                ')
    print('     AAAAA    MMMM      MMMM    AAAAA    LLL     GGGGGGGGG    AAAAA    MMMM      MMMM                ')
    print('    AAA AAA   MMMMM    MMMMM   AAA AAA   LLL     GGG   GGG   AAA AAA   MMMMM    MMMMM                ')
    print('   AAA   AAA  MMMMMM  MMMMMM  AAA   AAA  LLL     GGG   GGG  AAA   AAA  MMMMMM  MMMMMM                ')
    print('  AAA     AAA MMM MMMMMM MMM AAA     AAA LLL     GGG   GGG AAA     AAA MMM MMMMMM MMM     /^ ^\      ')
    print('  AAAAAAAAAAA MMM  MMMM  MMM AAAAAAAAAAA LLL      GGGGGGGG AAAAAAAAAAA MMM  MMMM  MMM    / 0 0 \     ')
    print('  AAA     AAA MMM   MM   MMM AAA     AAA LLL           GGG AAA     AAA MMM   MM   MMM    V\ Y /V     ')
    print('  AAA     AAA MMM        MMM AAA     AAA LLL           GGG AAA     AAA MMM        MMM    / - \       ')
    print('  AAA     AAA MMM        MMM AAA     AAA LLLLLLL  GGGGGGGG AAA     AAA MMM        MMM   /     |      ')
    print('  AAA     AAA MMM        MMM AAA     AAA LLLLLLL GGGGGGGGG AAA     AAA MMM        MMM   V__) ||      ')
    print('  -----------------------------------------------------------------------------------                ')
    print('  © Jasper A. Vrugt, University of California Irvine & GPT-4 OpenAI''s language model')
    print('    _________________________________________________________________________')
    print('    Version 2.0, Feb. 2025, Beta-release: MATLAB implementation is benchmark ')
    print("\n")

    return AMALGAMPar, Par_info, options, T_start


def AMALGAM_calc_setup(AMALGAMPar, fname, options, plugin):
    # ####################################################################### #
    # Sets up sequential / parallel  computational environment                #
    #  SYNOPSIS                                                               #
    #   [AMALGAMPar,f_handle] = AMALGAM_calc_setup(AMALGAMPar,fname, ...      #
    #       options,plugin)                                                   #
    # ####################################################################### #
    
    base_dir = None

    # Set up parallel execution based on options
    if options['parallel'] == 'no':
        AMALGAMPar['CPU'] = 1  # Use 1 CPU (processor)
    elif options['parallel'] == 'yes':
        # How many available workers? 
        workers = mp.cpu_count() 
        if workers > AMALGAMPar['N']:
            workers = AMALGAMPar['N']
        
        AMALGAMPar['CPU'] = workers
        
        # Write parallelization status to file
        with open('warning_file.txt', 'a+') as fid:
            msg = (f'AMALGAM PARALLEL: Pool opened with {AMALGAMPar["CPU"]} workers for a population of {AMALGAMPar["N"]} individuals\n')
            print(msg)
            fid.write(msg)
        
        # Handle I/O directories for parallel execution
        if options['IO'] == 'yes':
            if os.name in ['nt', 'posix']:  # For Windows or Unix systems

                example_dir = os.getcwd()                                   # Current working directory
                base_dir = os.path.join(example_dir, "worker_dirs")         # Base directory for worker directories
                # copy all files to model files
                temp_dir = os.path.join(os.getcwd(), 'model_files')         # Copy files to temporary directory
                if not os.path.exists(temp_dir):
                    os.makedirs(temp_dir)

                # Copy files to temporary directory
                for filename in os.listdir(example_dir):
                    file_path = os.path.join(example_dir, filename)
                    if os.path.isfile(file_path):  # Ensure only files are copied
                        shutil.copy(file_path, temp_dir)

                model_files_dir = os.path.join(example_dir, "model_files")  # Directory with all model files

            # Ensure the model files directory exists
            if not os.path.exists(model_files_dir):
                print(f"Model files directory '{model_files_dir}' does not exist!")
                
                return

            # Step 1: Initial setup, copy model files into each worker directory (done only once)
            if not os.path.exists(base_dir):
                os.makedirs(base_dir)

            for task_id in range(AMALGAMPar['CPU']):
                worker_dir = os.path.join(base_dir, f'worker_{task_id}')
                # Only copy the model files to the worker directories if not already done
                if not os.path.exists(worker_dir):
                    copy_model_files(model_files_dir, worker_dir)

    ## Function handle = dynamic function [not pickable, but ok with multiprocessing]
    func_handle = get_function_handle(fname)

    return AMALGAMPar, func_handle, base_dir


def AMALGAM_initialize(AMALGAMPar, Par_info, plugin):
    # ####################################################################### #
    #                                                                         #
    #  SYNOPSIS                                                               #
    #   [AMALGAMPar,Par_info,X,p_rm,PS,Z,output] = AMALGAM_initialize( ...    #
    #       AMALGAMPar,Par_info )                                             #
    #  where                                                                  #
    #                                                                         #
    # ####################################################################### #
    
    AMALGAMPar['q'] = len(AMALGAMPar['rec_methods'])            # Number of recombination methods
    p_rm = (1 / AMALGAMPar['q']) * np.ones(AMALGAMPar['q'])     # Selection probabilities for recombination methods
    AMALGAMPar['p0'] = max(2 / AMALGAMPar['N'], 5e-2)           # Minimum probability for recombination methods
    
    output = PS = {}
    output['p_rm'] = np.full((AMALGAMPar['T']+1, AMALGAMPar['q'] + 1), np.nan)  # Initialize matrix for p_rm
    output['IGD'] = np.full((AMALGAMPar['T']+1, 2), np.nan)                     # Initialize matrix for inverse generational distance
    output['p_rm'][0, :AMALGAMPar['q'] + 1] = np.concatenate(([0], p_rm))       # Store p_rm for recombination methods
    
    # Initialize Particle Swarm Optimization (PSO) if used
    if 'ps' in AMALGAMPar['rec_methods']:
        PS['v'] = (1 / 5) * np.random.uniform(-1, 1, (AMALGAMPar['N'], AMALGAMPar['d'])) * (Par_info['max'] - Par_info['min'])

    # Initialize the population
    X = np.full((AMALGAMPar['N'], AMALGAMPar['d']), np.nan)
    
    # Generate the initial population based on the specified method
    if Par_info['initial'] == 'uniform':
        X = Par_info['min'] + np.random.rand(AMALGAMPar['N'], AMALGAMPar['d']) * (Par_info['max'] - Par_info['min'])
    elif Par_info['initial'] == 'latin':
        X = LH_sampling(Par_info['min'], Par_info['max'], AMALGAMPar['N'])
    elif Par_info['initial'] == 'normal':
        X = np.tile(Par_info['mu'], (AMALGAMPar['N'], 1)) + np.random.randn(AMALGAMPar['N'], AMALGAMPar['d']) @ np.linalg.cholesky(Par_info['cov'])
    elif Par_info['initial'] == 'prior':
        if Par_info['u'] == 'yes':  # Univariate prior distribution
            for qq in range(AMALGAMPar['d']):
                for zz in range(AMALGAMPar['N']):
                    X
        else:  # Multivariate prior distribution
            for zz in range(AMALGAMPar['N']):
                X
    elif Par_info['initial'] == 'user':
        for zz in range(AMALGAMPar['N']):
            X[zz, :] = Par_info['x0'][zz, :]
    else:
        raise ValueError("AMALGAM_initialize:Unknown initial sampling method")

    # Boundary handling
    if 'boundhandling' in Par_info:
        X, v = Boundary_handling(X, Par_info)
    else:
        v = np.ones(AMALGAMPar['N'], dtype=bool)  # Initialize v for boundary checking (not used in the original)

    # BMA model training if applicable
    if 'unit_simplex' in Par_info:
        wght_sum = np.sum(X[:AMALGAMPar['N'], :int(plugin['BMA']['K'])], axis=1)
        X[:, :int(plugin['BMA']['K'])] = X[:, :int(plugin['BMA']['K'])] / wght_sum[:, np.newaxis]  # Normalize weights in the unit simplex

    # Transform to discrete space if required
    if 'steps' in Par_info:
        X = Discrete_space(X, Par_info)

    # Initialize the result matrix Z
    Z = np.full(((1 + AMALGAMPar['T']) * AMALGAMPar['N'] // AMALGAMPar['K'], AMALGAMPar['d'] + AMALGAMPar['m']), np.nan)
    
    return AMALGAMPar, Par_info, X, p_rm, PS, Z, output


def AMALGAM_rank(FQ, options=None):
    """
    Pareto ranking and crowding distance computation of solutions matrix FQ.

    Args:
    FQ (numpy.ndarray): 2D array of objective function values (size: M x m), 
                         where M is the number of solutions and m is the number of objectives.
    options (dict, optional): Algorithmic settings for density and ranking methods.
                              Default is None, in which case default settings will be used.

    Returns:
    RQ (numpy.ndarray): 1D array with the rank of each solution in FQ.
    dQ (numpy.ndarray): 1D array with the crowding distance or strength Pareto values.
    """
    if options is None:
#        options = {'density': 'crowding', 'ranking': 'matlab'}
        options = {'density': 'crowding'}

    M, m = FQ.shape  # Number of solutions and objective functions
    RQ = np.nan * np.ones(M)  # Initialize vector for ranks
    FQ_min = np.min(FQ, axis=0)  # Minimum objective function values

    # Initialize Pareto fronts and dominance structures
    Fr = [[]]  # Pareto-optimal fronts
    Sp = [[] for _ in range(M)]  # Points a particular point dominates
    nP = np.zeros(M, dtype=int)  # Number of points which dominate each point

    # Loop over all elements to determine the dominance
    for p in range(M):
        # Boolean indexing: check which points p dominates (less than or equal in all objectives)
        dominates = np.all(FQ[p, :] <= FQ, axis=1)  # Boolean array where p dominates other points
        strictly_dominates = np.any(FQ[p, :] < FQ, axis=1)  # Strict domination condition
        id = np.where(dominates & strictly_dominates)[0]  # Points dominated by p

        if len(id) > 0:
            Sp[p] = id
        # Find which points dominate p (strict domination condition)
        dominates_p = np.all(FQ <= FQ[p, :], axis=1)
        strictly_dominates_p = np.any(FQ < FQ[p, :], axis=1)
        id = np.where(dominates_p & strictly_dominates_p)[0]
        nP[p] = len(id)   # original translation GPT
        # nP[p] = nP[p] + len(id)
        if nP[p] == 0:  # p is in the first front
            Fr[0].append(p)

    # Determine the Pareto rank
    i = 0
    while len(Fr[i]) > 0:
        NextFr = []  # Next front
        for p in Fr[i]:
            q = Sp[p]
            # Ensure that q is an array of indices and that nP[q] can be accessed properly
            nP[q] -= 1  # Decrease the dominance count for each point in the front
            # Fix the indexing to avoid the boolean mismatch
            indices = np.where(nP[q] == 0)[0]  # Get the indices where nP[q] == 0
            if len(indices) != 0:           # Added this line
                NextFr.extend(q[indices])   # Add the new front members
        i += 1
        Fr.append(NextFr)

    # Assign ranks based on the front structure
    for j in range(i):
        RQ[Fr[j]] = j + 1  # Assign ranks to solutions

    # Calculate crowding distance or strength Pareto based on the density method
    if options['density'] == 'crowding':
        dQ = np.nan * np.ones(M)
        for j in range(int(np.max(RQ))):
            id = np.where(RQ == j + 1)[0]
            R_s = FQ[id, :]
            N_s = len(id)
            Cr_d = np.zeros(N_s)
            for i in range(m):
                sorted_indices = np.argsort(R_s[:, i])
                Cr_d[sorted_indices[0]] = np.inf
                Cr_d[sorted_indices[-1]] = np.inf
                for z in range(1, N_s - 1):
                    Cr_d[sorted_indices[z]] += (R_s[sorted_indices[z + 1], i] - R_s[sorted_indices[z - 1], i])
            dQ[id] = Cr_d
    elif options['density'] == 'strength':
        n_dom = np.zeros(M)
        for j in range(M):
            n_dom[j] = np.sum(np.all(FQ[j, :] < FQ, axis=1))
        dQ = M / n_dom  # Strength Pareto
    else:
        raise ValueError("Unknown density estimation method")

    return RQ, dQ, FQ_min


def AMALGAM_restart(file_name, Func_name, Ftrue):
    """
    Restart function to complete the desired number of generations.

    Parameters:
    fname : str
        The file name containing the saved AMALGAM data to restart from.

    Returns:
    tuple : (AMALGAMPar, Par_info, options, PS, X, Z, FX, YX, output, FX_min, RX, dX, T_start)
        The updated AMALGAM structure, function handle, data arrays, and other variables.
    """
    # Try-except block to handle loading failure
    try:
        # Load the restart data from the file `file_name`
        data = np.load(file_name + '.npy', allow_pickle=True)           # with shelve.open(file_name, 'r') as file:
        loaded_data = data.item()
        AMALGAMPar = loaded_data['AMALGAMPar']                          #   AMALGAMPar = file['AMALGAMPar']
        Par_info = loaded_data['Par_info']                              #   Par_info = file['Par_info']
        options = loaded_data['options']                                #   options = file['options']
        PS = loaded_data['PS']                                          #   PS = file['PS']
        output = loaded_data['output']                                  #   output = file['output']
        X = loaded_data['X']                                            #   X = file['X']
        Z = loaded_data['Z']                                            #   Z = file['Z']
        FX = loaded_data['FX']                                          #   FX = file['FX']
        YX = loaded_data['YX']                                          #   YX = file['YX']
        FX_min = loaded_data['FX_min']                                  #   FX_min = file['FX_min']
        RX = loaded_data['RX']                                          #   RX = file['RX']
        dX = loaded_data['dX']                                          #   dX = file['dX']
        p_rm = loaded_data['p_rm']                                      #   p_rm = file['p_rm']
        t = loaded_data['t']                                            #   t = file['t']
        ct = loaded_data['ct']                                          #   ct = file['ct']
        
        # Load the content from the pickle file
        #  with open("{}.pkl".format(file_name), 'rb') as file:
        #      plugin = pickle.load(file)

    except FileNotFoundError:
        # Handle the case when the file does not exist
        error_message = f"AMALGAM_PACKAGE ERROR: Cannot restart --> File {fname} does not exist. Next run, to avoid this problem, set field 'save' of structure options to 'yes'"
        raise FileNotFoundError(error_message)

    # Open warning file
    with open('warning_file.txt', 'a+') as fid:
        # Check whether the previous run was aborted early or not
        if t < AMALGAMPar['T']:
            # Finish previous budget
            warning_message = f"AMALGAM RESTART: Starting with t = {t} but still using old budget of AMALGAMPar.T = {AMALGAMPar['T']}\n"
            print(warning_message)
            fid.write(warning_message)
        else:
            # Assign new budget and run
            fid.write('-------------- AMALGAM warning file --------------\n')
            t_file = 'T.txt'
            
            if os.path.exists(t_file):
                # If T.txt exists, check its content
                with open(t_file, 'r') as f:
                    T_new = checkfile_T_AMALGAM(t_file)
                
                warning_message = f"AMALGAM RESTART: User has requested/listed {T_new} additional generations in file 'T.txt'\n"
                print(warning_message)
                fid.write(warning_message)

                final_message = f"AMALGAM RESTART: Initial t = {t} and completing {T_new} additional generations so that AMALGAMPar.T = {AMALGAMPar['T'] + T_new}\n"
                print(final_message)
                fid.write(final_message)

            else:
                # If T.txt does not exist, double the previous budget
                warning_message = "AMALGAM RESTART: Did not locate file 'T.txt' in respective example directory so per default we double the 'previous' budget of generations\n"
                print(warning_message)
                fid.write(warning_message)
                
                T_new = AMALGAMPar['T']
                final_message = f"AMALGAM RESTART: Initial t = {t} and completing {T_new} additional generations so that AMALGAMPar.T = {AMALGAMPar['T'] + T_new}\n"
                print(final_message)
                fid.write(final_message)

            # Add T_new * AMALGAMPar['N'] rows to Z with nan values for T_new additional generations
            Z = np.pad(Z, ((0, T_new * AMALGAMPar['N']), (0, 0)), mode='constant', constant_values=np.nan)
            # Add T_new lines to p_rm
            output['p_rm'] = np.pad(output['p_rm'], ((0, T_new), (0, 0)), mode='constant', constant_values=np.nan)
            if len(Ftrue) > 0:
                # Add T_new lines to IGD
                output['IGD'] = np.pad(output['IGD'], ((0, T_new), (0, 0)), mode='constant', constant_values=np.nan)
            AMALGAMPar['T'] += T_new

    # Define starting value of T
    T_start = t + 1

    # Setup parallel computing framework or not (this part would need further context to fully translate)
    AMALGAMPar, func_handle, base_dir = AMALGAM_calc_setup(AMALGAMPar, Func_name, options, plugin)

    return AMALGAMPar, Par_info, func_handle, options, PS, X, Z, FX, YX, p_rm, output, FX_min, RX, dX, ct, base_dir, T_start


def AMALGAM_distribution(AMALGAMPar, p_rm):
    """
    Determine load distribution among recombination methods.
    
    Parameters:
    AMALGAMPar : dict
        A dictionary containing the parameters, including:
        - 'q': the number of recombination methods
        - 'N': the total number of points
    p_rm : list or numpy array
        A list or array representing the probability distribution of recombination methods.
        
    Returns:
    id : numpy array
        An array indicating which recombination method is responsible for each point.
    id_rm : list of numpy arrays
        A list of arrays, where each element is a list of indices corresponding to a recombination method.
    """
    # Ensure that p_rm is a horizontal vector (1D array)
    p_rm = np.asarray(p_rm).flatten()
    
    # Randomly sample recombination methods for each point based on probabilities p_rm
    id = np.random.choice(np.arange(1, AMALGAMPar['q'] + 1), size=AMALGAMPar['N'], p=p_rm)
    
    # Create a list of indices for each recombination method
    id_rm = [np.where(id == j)[0] for j in range(1, AMALGAMPar['q'] + 1)]
    
    return id, id_rm


def AMALGAM_population(AMALGAMPar, options, X, G, FX, FG, id):
    """
    Selects the new population based on current offspring and parents.
    
    Parameters:
    AMALGAMPar : dict
        A dictionary containing the parameters, including:
        - 'N': the population size
        - 'd': the dimensionality of the problem
        - 'm': the number of objectives
    options : dict
        A dictionary of options for the algorithm, including selection details.
    X : numpy array
        The parent population.
    G : numpy array
        The offspring population.
    FX : numpy array
        The objective function values of the parent population.
    FG : numpy array
        The objective function values of the offspring population.
    id : numpy array
        An array of recombination method indices.
        
    Returns:
    Xn : numpy array
        The new population (after selection).
    FXn : numpy array
        The objective function values of the new population.
    RXn : numpy array
        The ranks of the new population.
    dXn : numpy array
        The crowding distances of the new population.
    id_N : numpy array
        The indices of the new population.
    id : numpy array
        The recombination method indices for the new population.
    """
    
    # Combine parent and offspring populations
    Q = np.vstack([X, G])
    FQ = np.vstack([FX, FG])
    
    # Rank and calculate crowding distances
    RQ, dQ, _ = AMALGAM_rank(FQ, options)  # AMALGAM_rank needs to be defined elsewhere
    
    # Indices of recombination methods
    I_alg = np.hstack([np.zeros(AMALGAMPar['N']), id])
    
    # Selection based on rank
    RQ_max = int(np.max(RQ))  # Maximum rank of Q
    n_rnk = np.full((RQ_max, 1), np.nan)

    for r in range(1, RQ_max + 1):
        id_r = np.where(RQ == r)[0]
        n_rnk[r - 1] = len(id_r)  # Number of points with rank "r"
    
    tot_rnk = np.cumsum(n_rnk)  # Cumulative sum of the number of points with each rank
    idx = np.where(tot_rnk <= AMALGAMPar['N'])[0]  # Find ranks that fit within the population size

    if idx.size > 0:
        r_max = idx[-1]  # Maximum rank
        id_R = np.where(RQ <= r_max+1)[0]  # Select all points with rank <= r_max
        #n_lft = AMALGAMPar['N'] - tot_rnk[r_max - 1]  # Remaining spots to fill
        n_lft = int(AMALGAMPar['N'] - tot_rnk[r_max])  # Remaining spots to fill
        id_lft = np.where(RQ == r_max + 2)[0]  # Points with the next rank
        # Sort based on crowding distance (descending)
        sorted_lft = np.argsort(-dQ[id_lft])
        id_C = id_lft[sorted_lft[:n_lft]]  # Select points with the largest crowding distances
        # Combine the selected points
        id_N = np.hstack([id_R, id_C])
    else:
        id_R = np.where(RQ == 1)[0]  # Select only points with rank 1
        sorted_R = np.argsort(-dQ[id_R])  # Sort based on crowding distance (descending)
        id_N = id_R[sorted_R[:AMALGAMPar['N']]]  # Select the top N points based on crowding distance
    
    # Extract the new population
    Xn = Q[id_N, :AMALGAMPar['d']]  # New population
    FXn = FQ[id_N, :AMALGAMPar['m']]  # Objective function values of the new population
    RXn = RQ[id_N]  # Ranks of the new population
    dXn = dQ[id_N]  # Crowding distances of the new population
    id = I_alg[id_N]  # Indices of the recombination methods for the new population
    
    return Xn, FXn, RXn, dXn, id_N, id


def AMALGAM_load(AMALGAMPar, p_rm, id):
    """
    Determines the new number of points to generate with individual algorithms.
    
    Parameters:
    AMALGAMPar : dict
        A dictionary containing the parameters, including:
        - 'q': the number of recombination methods
        - 'p0': the minimum selection probability
    p_rm : numpy array
        The initial probabilities for each recombination method.
    id : numpy array
        The recombination method indices for the new population.
        
    Returns:
    p_rm : numpy array
        The updated probabilities for each recombination method.
    """
    
    # Initialize an array to store the number of points for each recombination method
    n_rm = np.full(AMALGAMPar['q'], np.nan)  # Initialize with NaNs

    # Calculate the total number of points made it to the new population (id > 0)
    T = np.sum(id > 0)  # Number of points that made it to the new population
    
    # For each recombination method, calculate the proportion of points from that method
    for j in range(1, AMALGAMPar['q']+1):
        n_rm[j-1] = np.sum(id == j) / T  # Proportion of points from recombination method j
    
    # Scale the contribution of each method by its selection probability
    p = n_rm / p_rm  # Scale by current probabilities
    
    # Normalize the probabilities and ensure the minimum selection probability (AMALGAMPar.p0)
    pn = np.maximum(p / np.sum(p), AMALGAMPar['p0'])
    
    # Normalize again to ensure the probabilities sum to 1
    p_rm = pn / np.sum(pn)

    return p_rm


def AMALGAM_children(AMALGAMPar, Par_info, X, FX, RX, dX, PS, id_rm, plugin):
    """
    This function uses different recombination methods to generate children.

    Parameters:
    AMALGAMPar : dict
        Contains parameters such as the number of children and the dimension.
    Par_info : dict
        Contains additional parameter info, including boundaries, unit simplex, etc.
    X : numpy.ndarray
        Current population of solutions.
    FX : numpy.ndarray
        Objective function values corresponding to the population.
    RX : numpy.ndarray
        Additional information needed for some methods (e.g., random numbers).
    dX : numpy.ndarray
        Difference matrix for use in recombination methods.
    PS : dict
        Information related to Particle Swarm optimization.
    id_rm : list
        List of indices for individuals to be removed or used for creating children.
        
    Returns:
    G : numpy.ndarray
        The generated children.
    """
    
    # Initialize offspring as NaN
    G = np.full((AMALGAMPar['N'], AMALGAMPar['d']), np.nan)

    # <><><><><><><><><><><> NSGA-II GENETIC ALGORITHM <><><><><><><><><><><><>
    if 'ga' in AMALGAMPar['rec_methods']:
        run_GA = AMALGAMPar['rec_methods'].index('ga')
        id_GA = id_rm[run_GA]  # Indices of population
        G[id_GA, :] = NSGA(AMALGAMPar, Par_info, X, FX, id_GA, dX)  # Apply NSGA

    # <><><><><><><><><><> ADAPTIVE METROPOLIS ALGORITHM <><><><><><><><><><><>
    if 'am' in AMALGAMPar['rec_methods']:
        run_AM = AMALGAMPar['rec_methods'].index('am')
        id_AM = id_rm[run_AM]
        n_AM = len(id_AM)  # Number of children AM must create
        G[id_AM, :] = AM(AMALGAMPar, X, RX, n_AM, id_AM)  # Apply Adaptive Metropolis
    
    # <><><><><><><><><><><><> PARTICLE SWARM METHOD <><><><><><><><><><><><><>
    if 'ps' in AMALGAMPar['rec_methods']:
        run_PS = AMALGAMPar['rec_methods'].index('ps')
        id_PS = id_rm[run_PS]
        n_PS = len(id_PS)  # Number of children PS must create
        G[id_PS, :], PS = PSO(AMALGAMPar, X, PS, n_PS, id_PS)  # Apply Particle Swarm
#        PS['v'] = V  # Update velocity

    # <><><><><><><><><><><><> DIFFERENTIAL EVOLUTION <><><><><><><><><><><><><
    if 'de' in AMALGAMPar['rec_methods']:
        run_DE = AMALGAMPar['rec_methods'].index('de')
        id_DE = id_rm[run_DE]
        n_DE = len(id_DE)  # Number of children DE must create
        G[id_DE, :] = DE(AMALGAMPar, X, n_DE, id_DE)  # Apply Differential Evolution
    
    # Boundary handling
    if 'boundhandling' in Par_info:
        G, v = Boundary_handling(G, Par_info)  # Adjust according to boundary handling
    else:
        v = np.ones(AMALGAMPar['N'], dtype=bool)  # Default boundary handling
    
    # BMA model training (if unit simplex is specified)
    if 'unit_simplex' in Par_info:
        w_sum = np.sum(G[:, :int(plugin['BMA']['K'])], axis=1)  # Sum of BMA weights
        G[:, :int(plugin['BMA']['K'])] = G[:, :int(plugin['BMA']['K'])] / w_sum[:, None]  # Normalize weights
    
    # Transform to discrete space (if steps are defined)
    if 'steps' in Par_info:
        G = Discrete_space(G, Par_info)  # Transform to discrete space
    
    return G, PS


def AMALGAM_calc_FX(X, AMALGAMPar, options, func_handle, base_dir, plugin, printed_warnings, verbose = 0):
    """
    Evaluate user-supplied function and return objective function values or model simulations if so desired.

    Parameters:
    X : np.ndarray
        The parameter vectors (size N x d).
    AMALGAMPar : dict
        Dictionary containing algorithm parameters such as 'm' (number of objectives), 'CPU', etc.
    options : dict
        A dictionary containing options such as 'modout' and 'IO'.
    verbose : int, optional
        If set to 1, prints progress.

    Returns:
    FX : np.ndarray
        The objective function values for each parameter vector.
    Y : np.ndarray
        Model simulations (if 'modout' is 'yes'), otherwise an empty array.
    """
    N, d = X.shape                  # Number of parameter vectors and dimensions
    m = AMALGAMPar['m']             # Number of objective functions
    FX = np.full((N, m), np.nan)    # Preallocate objective function values
    Y = None                        # Preallocate model simulation output

    # Sequential evaluation - CORRECT in Python
    if AMALGAMPar['CPU'] == 1:
        for ii in range(N):
            # If function returns more than 1 output, then it will go to *rest, but be careful
            if plugin is None:
                results = func_handle(X[ii, :])
            else:
                results = func_handle(X[ii, :], plugin)
            # results can consists or more than one array depending on # return arguments Func_name
            if isinstance(results, (tuple, list)):
                FX[ii, :] = results[0]
                # Now check whether we have to store model output
                if options['modout'] == 'yes':
                    if ii == 0:                                     # Initialize Y on first iteration
                        # Z = np.array([x[1] for x in results])     # Extract 2nd element
                        Z = results[1]                              # Extract 2nd element
                        Y = np.full((N, len(Z)), np.nan)            # Initialize Y 
                        Y[0, :] = Z                                 # Populate first row of Y
                    else:
                        Y[ii, :] = results[1]                       # Subsequent iterations
                else:
                    Y = None                                        # Y is returned but modout is 'no'
            else:
                FX[ii,] = results
                Y = None
            if verbose:
                # Print progress if verbose flag is set
                prt_progress(AMALGAMPar, N)        

    # Parallel evaluation - CORRECT in Python if IO = Yes and IO = No
    elif AMALGAMPar['CPU'] > 1:
        task_ranges = distribute_tasks(N, AMALGAMPar['CPU'])      # Task distribution for each worker (divide work)
        # task_ranges = [(i * (N // DREAMPar['CPU']), (i + 1) * (N // DREAMPar['CPU'])) for i in range(DREAMPar['CPU'])]
        if isinstance(plugin, dict):                            # Convert to regular numpy array so that it can be shared with workers [= pickable]
            plugin = convert_memoryview_to_array(plugin)
            
        with mp.Pool(processes=AMALGAMPar['CPU']) as pool:
            results = pool.starmap(worker_task, [(worker_id, start_idx, end_idx, X, func_handle, plugin, base_dir)
                                                    for worker_id, (start_idx, end_idx) in enumerate(task_ranges)])

        if isinstance(results, (tuple, list)):      
            # Unpack the results from the workers
            ct = 0
            i = 0
            for worker_result in results:
                # Determine how many variables [objects] were received from function handle
                nvar = sum(isinstance(item, np.ndarray) for sublist in worker_result for item in sublist)
                if N > AMALGAMPar['CPU']:     ## must divide by number of X's evaluated by each worker
                    nvar = nvar/(task_ranges[i][1] - task_ranges[i][0])
                i = i + 1
                if isinstance(worker_result, list):         ## This is when N > DREAMPar['CPU']
                    if nvar == 2:    ## 2 outputs, FX and Y    
                        for fx, y in worker_result:
                            FX[ct, :] = fx.flatten() 
                            if options['modout'] == 'no':
                                warning_msg = f"AMALGAM_calc_FX WARNING: Did not expect two output arguments from the function. Setting Y to None as options['modout'] == 'no'."
                                if warning_msg not in printed_warnings:
                                    print(warning_msg)
                                    printed_warnings.add(warning_msg)
                            else:
                                if ct == 0:
                                    Y = np.full((N,len(y)),np.nan)
                                Y[ct, :] = y.flatten()
                            # update counter                                
                            ct = ct + 1
    
                    elif nvar == 1:  ## 1 output, FX
                        for fx in worker_result:
                            FX[ct, :] = fx.flatten()
                            ct = ct + 1
                            if options['modout'] == 'yes':
                                warning_msg = f"AMALGAM_calc_FX WARNING: Did not expect one output argument from the function as options['modout'] == 'yes'. Setting Y to None"
                                if warning_msg not in printed_warnings:
                                    print(warning_msg)
                                    printed_warnings.add(warning_msg)
    
                    else:           ## 0 or more than 2 outputs
                        warning_msg = f"AMALGAM_calc_FX ERROR: No output arguments from the function. Setting FX and Y to None."
                        if warning_msg not in printed_warnings:
                            print(warning_msg)
                            printed_warnings.add(warning_msg)
 
                elif isinstance(worker_result, tuple):      ## This is when N = DREAMPar['CPU']
                    if nvar == 2:    ## 2 outputs, FX and Y
                        fx, y = worker_result
                        if options['modout'] == 'no':
                            warning_msg = f"AMALGAM_calc_FX WARNING: Did not expect two output arguments from the function. Setting Y to None as options['modout'] == 'no'."
                            y = []
                            if warning_msg not in printed_warnings:
                                print(warning_msg)
                                printed_warnings.add(warning_msg)
 
                    elif nvar == 1:  ## 1 output, FX
                        fx = worker_result
                        y = []
                        if options['modout'] == 'yes':
                            warning_msg = f"AMALGAM_calc_FX WARNING: Did not expect one output argument from the function as options['modout'] == 'yes'. Setting Y to None"
                            if warning_msg not in printed_warnings:
                                print(warning_msg)
                                printed_warnings.add(warning_msg)

                    else:           ## 0 or more than 2 outputs
                        warning_msg = f"AMALGAM_calc_FX ERROR: No output arguments from the function. Setting FX and Y to None."
                        if warning_msg not in printed_warnings:
                            print(warning_msg)
                            printed_warnings.add(warning_msg)

                    if options['modout'] == 'yes':
                        FX[ct, :] = fx.flatten()
                        if ct == 0:
                            Y = np.full((N,len(y)),np.nan)
                        Y[ct, :] = y.flatten()
                    else:
                        FX[ct, :] = fx.flatten()
                    # update counter        
                    ct = ct + 1  

        # If there is only one return value, unpack it directly
        elif not isinstance(results, (tuple, list)):
            FX = results
            Y = None
    
    if verbose:
        print("\nModel simulation ... done")

    return FX, Y


def AMALGAM_end(AMALGAMPar, options, base_dir):
    """
    Finalizes the AMALGAM process, closes workers, deletes worker directories, 
    and writes the final warning message to a file.
    
    Parameters:
    AMALGAMPar : dict
        Contains the parameters, specifically the CPU count.
    options : dict
        Contains options such as IO operations.
    """
    # Close workers if there are multiple CPUs
    if AMALGAMPar['CPU'] > 1:
        # Close parallel pool: check in Python
        pass

        if options['IO'] == 'yes':  # If IO writing is enabled, remove directories
            # Step 3: Clean up (delete) the worker directories after all generations are done
            cleanup_worker_directories(base_dir, AMALGAMPar['CPU'])

    # Open the warning_file.txt file in append mode
    with open('warning_file.txt', 'a+') as fid:
        # Write final line of warning file
        fid.write('----------- End of AMALGAM warning file ----------\n')

    # Optionally open the warning file in the default editor for the platform
    # if platform.system() in ['Windows', 'Darwin']:  # Windows or macOS
    #     os.system('notepad warning_file.txt' if platform.system() == 'Windows' else 'open warning_file.txt')
    return


def AMALGAM_postproc(AMALGAMPar, Par_info, options, X, FX, output, Fpareto=None, YX = None, Z = None):
    # Default behavior for YX
#    if YX is None:
#        YX = []

    # Print wait statement to the screen
    print('AMALGAM PLOTTING: PLEASE WAIT ...')
    
    # Define name of program
    n_program = 'AMALGAM'
    
    # Define name of figures file
    file_name = f'{n_program}_figures.pdf'
    
    # Determine screen size (using matplotlib to get screen dimensions)
    monitor = get_monitors()[0]
    screen_width = monitor.width
    screen_height = monitor.height
    x_mult = screen_width / 1920
    y_mult = screen_height / 1080
    t_mult = min(x_mult, y_mult)

    # Define fontsize for figures
    fontsize_xylabel = 16 * t_mult
    fontsize_axis = 16 * t_mult
    fontsize_legend = 14 * t_mult
    fontsize_text = 14 * t_mult
    fontsize_title = 18 * t_mult
    markersize_symbols = 5
    fontsize_titlepage = 20 * t_mult
    markersize_legend = 10 * t_mult

    symbol = ['r', 'b', 'g', 'c', 'k', 'm', 'y']
    # prepare parameter names
    maxbins = 25; str_par = []

    if 'names' in Par_info:
        for name in Par_info['names']:
            str_par.append(f"$\\;{name}$")
            # Now create parameter strings for the Tables
        str_table = [s.replace('\text{', '').replace('}', '').replace('{', '').replace('\\', '') for s in Par_info['names']]
    else:
        for i in range(AMALGAMPar['d']):
            str_par.append(f"$\\;x_{{{i+1}}}$")
            # Prepare table string (for output formatting)
        str_table = [s.replace('$', '').replace('\\;', '') for s in str_par]

    max_length = max(len(name) for name in str_table)
    # Extract min and max values from Par_info
    Par_info['min'] = Par_info['min'][0, :]
    Par_info['max'] = Par_info['max'][0, :]

    # Rank the final population and focus on rank 1 solutions
    R, _, _ = AMALGAM_rank(FX, options)  # Assuming AMALGAM_rank is defined elsewhere
    idR = (R == 1)

    # Isolate parameter and objective function values of rank 1 solutions
    X1 = X[idR, :]
    FX1 = FX[idR, :]

    # Calculate mean, standard deviation, and median of Pareto solutions
    MEAN = np.mean(X1, axis=0)
    STD = np.std(X1, axis=0)
    MED = np.median(X1, axis=0)

    # Determine optimum of each objective function
    MAP = np.full((AMALGAMPar['m'], AMALGAMPar['d']), np.nan)
    for ii in range(AMALGAMPar['m']):
        sorted_indices = np.argsort(FX[:AMALGAMPar['N'], ii])
        MAP[ii, :] = X[sorted_indices[0], :]

    # Write the results to an output file
    with open('AMALGAM_output.txt', 'w') as fid:
        fid.write('--------------------------- AMALGAM output file ------------------ \n')
        fid.write('\n')
        fid.write(' Table 1: Parameter median, mean and standard deviation of rank 1 solutions \n')
        fid.write('          ================================================== \n')
        fid.write('             x             MEDIAN      MEAN      STD         \n')
        fid.write('          -------------------------------------------------- \n')
        # Print parameter statistics
        # fmt_1 = '            %-4s \t %7.3f    %7.3f    %7.3f\n'
        fmt_1 = f'          %-{max_length}s \t %7.3f    %7.3f    %7.3f\n'
        for j in range(AMALGAMPar['d']):
            fid.write(fmt_1 % (str_table[j], MED[j], MEAN[j], STD[j]))
        fid.write('          ================================================== \n')
        fid.write('\n')
        fid.write('\n')

        # Calculate the correlation matrix
        CORR = np.corrcoef(X1, rowvar=False)
        
        table_width = (max_length+3) * ( AMALGAMPar['d'] + 1)  # Each column is 7 characters wide
        top_line = '=' * table_width    # Line of '=' with the correct length
        # Print correlation matrix
        fid.write(' Table 2: Pearson''s correlation coefficients among parameters of rank 1 solutions \n')
        fid.write(f"          {top_line} \n")
        fid.write(f"{'':{max_length+17}}")
        for i in range(AMALGAMPar['d']):
            fid.write(f"{str_table[i]:{max_length+3}}")
        fid.write('\n')
        for i in range(AMALGAMPar['d']):
            fid.write(f"          {str_table[i]:{max_length+3}}")
            for j in range(AMALGAMPar['d']):
                fid.write(f"{CORR[i, j]:{max_length+3}.3f}")
            fid.write('\n')
        fid.write(f"          {top_line} \n")
        fid.write('\n')
        fid.write('------------------------- End AMALGAM output file ---------------- \n')

    with PdfPages(file_name) as pdf:

        ### Plot Empty Figure for PDF
        plt.figure(figsize=(12, 6))
        plt.plot([], [], 'ro')  # Empty plot
        plt.axis([0, 1, 0, 1])
        plt.gca().set_facecolor('w')
        plt.gca().set_xticklabels([])
        plt.gca().set_yticklabels([])        
        plt.gca().set_xticks([])
        plt.gca().set_yticks([])
        plt.text(0.3 * x_mult, 0.6 * y_mult, r'${\rm Visual \; results \; of \; AMALGAM \; toolbox}$', fontsize=fontsize_titlepage) #, ha='center', va='center')
        plt.text(0.27 * x_mult, 0.5 * y_mult, r'$\;\;\;\;\;{\rm Tables \; are \; not \; printed \; to \; PDF \; file}$', fontsize=fontsize_titlepage) #, ha='center', va='center') #, fontweight='bold')
        ax = plt.gca()  # Get current axis
        ax.spines['top'].set_visible(False)
        ax.spines['right'].set_visible(False)
        ax.spines['bottom'].set_visible(False)
        ax.spines['left'].set_visible(False)
        pdf.savefig()    
        plt.show()

        ### 3D Pareto Front Plot
        if AMALGAMPar['m'] == 3:
            fig = plt.figure(figsize=(12, 10))
            ax = fig.add_subplot(111, projection='3d')
            if len(Fpareto) > 0:
                ax.scatter(Fpareto[:, 0], Fpareto[:, 1], Fpareto[:, 2], color='k', s=50, label='True Pareto front')
            ax.scatter(FX1[:, 0], FX1[:, 1], FX1[:, 2], color='gray', s=50, label='Rank 1 Solution')
            ax.set_xlabel('$f_1$', fontsize=fontsize_xylabel)
            ax.set_ylabel('$f_2$', fontsize=fontsize_xylabel)
            ax.set_zlabel('$f_3$', fontsize=fontsize_xylabel)
            ax.set_title('AMALGAM: Three-dimensional plot of Pareto solution set in objective function space', fontsize=fontsize_title)
            # Add optimal points
            ax.scatter(FX1[np.argmin(FX1[:, 0]), 0], FX1[np.argmin(FX1[:, 0]), 1], FX1[np.argmin(FX1[:, 0]), 2], color=symbol[0], s=100, label='Optima 1')
            ax.scatter(FX1[np.argmin(FX1[:, 1]), 0], FX1[np.argmin(FX1[:, 1]), 1], FX1[np.argmin(FX1[:, 1]), 2], color=symbol[1], s=100, label='Optima 2')
            ax.scatter(FX1[np.argmin(FX1[:, 2]), 0], FX1[np.argmin(FX1[:, 2]), 1], FX1[np.argmin(FX1[:, 2]), 2], color=symbol[2], s=100, label='Optima 3')
            if len(Fpareto) > 0:
                legend_handles = [Line2D([0], [0], color='k', lw=0, marker='o', markersize=markersize_legend, label=f'True Pareto front')]
                legend_handles[1:] = [Line2D([0], [0], color='gray', lw=0, marker='o', markersize=markersize_legend, label=f'Rank 1 solution')]
                legend_handles[2:4] = [Line2D([0], [0], color=symbol[i], lw=0, marker='o', markersize=markersize_legend, label=r'$F_{{{}, \text{{opt}}}}$'.format(i+1)) for i in range(0,3)]
                symbol_ext = ['k','gray'] + symbol
            else:
                legend_handles = [Line2D([0], [0], color='gray', lw=0, marker='o', markersize=markersize_legend, label=f'Rank 1 solution')]
                legend_handles[1:3] = [Line2D([0], [0], color=symbol[i], lw=0, marker='o', markersize=markersize_legend, label=r'$F_{{{}, \text{{opt}}}}$'.format(i+1)) for i in range(0,3)]
                symbol_ext = ['gray'] + symbol 
                
            legend = ax.legend(handles=legend_handles, fontsize=fontsize_legend, loc='upper right')
            # Match text color to line color
            for z, text in enumerate(legend.get_texts()):
                text.set_color(symbol_ext[z])
  #          ax.legend(fontsize=fontsize_legend, loc='upper right')
            pdf.savefig()        
            plt.show()

        ### 2D Pareto Front Snapshots
        C = list(combinations(range(AMALGAMPar['m']), 2))  # All 2D combinations of objectives
        n_C = len(C)
        if n_C == 1:
            row = 1; col = 2
        else:
            row = 2; col = 3
        counter = 0
        for i in range(n_C):
            if counter == 0 or counter == row*col:
                fig, axs = plt.subplots(row, col, figsize=(15, 10))
                axs = axs.flatten()  # Flatten the 2D array of subplots to 1D for easier indexing
                counter = 0
            ax = axs[counter]
            d1, d2 = C[i]
            if len(Fpareto) > 0:
                ax.scatter(Fpareto[:, d1], Fpareto[:, d2], color='k', s=50, label='True Pareto front')
            ax.scatter(FX1[:, d1], FX1[:, d2], color='gray', s=50, label='Rank 1 Solution')
            # Plot optimal points
            ax.scatter(FX1[np.argmin(FX1[:, d1]), d1], FX1[np.argmin(FX1[:, d1]), d2], color=symbol[d1], s=100, label=r'$F_{{{}, \text{{opt}}}}$'.format(d1+1))
            ax.scatter(FX1[np.argmin(FX1[:, d2]), d1], FX1[np.argmin(FX1[:, d2]), d2], color=symbol[d2], s=100, label=r'$F_{{{}, \text{{opt}}}}$'.format(d2+1))
            ax.set_xlabel(f'$f_{{{d1+1}}}$', fontsize=fontsize_xylabel)
            ax.set_ylabel(f'$f_{{{d2+1}}}$', fontsize=fontsize_xylabel)
            d1d2 = np.array([d1,d2])
#            ax.legend(fontsize=fontsize_legend)
            if len(Fpareto) > 0:
                legend_handles = [Line2D([0], [0], color='k', lw=3, label=f'True Pareto front')]
                legend_handles[1:] = [Line2D([0], [0], color='gray', lw=0, marker='o', markersize=markersize_legend, label=f'Rank 1 solution')]
                legend_handles[2:3] = [Line2D([0], [0], color=symbol[d1d2[i]], lw=0, marker='o', markersize=markersize_legend, label=r'$F_{{{}, \text{{opt}}}}$'.format(d1d2[i]+1)) for i in range(0,2)]
                symbol_ext = ['k','gray'] + [symbol[key] for key in d1d2]
            else:
                legend_handles = [Line2D([0], [0], color='gray', lw=0, marker='o', markersize=markersize_legend, label=f'Rank 1 solution')]
                legend_handles[1:AMALGAMPar['m']] = [Line2D([0], [0], color=symbol[d1d2[i]], lw=0, marker='o', markersize=markersize_legend, label=r'$F_{{{}, \text{{opt}}}}$'.format(d1d2[i]+1)) for i in range(0,2)]
                symbol_ext = ['gray'] + [symbol[key] for key in d1d2]
                
            legend = ax.legend(handles=legend_handles, fontsize=fontsize_legend, loc='upper right')
            # Match text color to line color
            for z, text in enumerate(legend.get_texts()):
                text.set_color(symbol_ext[z])  
            if counter == 0:
                 # Set the title of the figure
                fig.suptitle('AMALGAM: Two-dimensional snapshots of Pareto samples in objective function space', fontsize=fontsize_title, y=0.92)
            # increment counter                
            counter += 1
            # check whether we should save/show figure or not
            if counter == row*col or i == n_C - 1:
                # delete empty subplots
                for ax in axs:
                    if not ax.has_data():  # If the axis has no data, remove it
                        fig.delaxes(ax)
                pdf.savefig()
                plt.show()

        ### Correlation Matrix Plot
        if AMALGAMPar['d'] > 1:
            fig, ax = plt.subplots(figsize=(8, 6))
            cax = ax.matshow(CORR, cmap='coolwarm')
            plt.colorbar(cax)
            ax.set_xticks(np.arange(AMALGAMPar['d']))
            ax.set_yticks(np.arange(AMALGAMPar['d']))
            ax.set_xticklabels(str_par)
            ax.set_yticklabels(str_par)
            plt.title('AMALGAM: Map of correlation coefficients of Pareto parameter samples', fontsize=fontsize_title*3/4)
            pdf.savefig()        
            plt.show()
        elif AMALGAMPar['d'] == 1:
            print("AMALGAM WARNING: Cannot plot map with bivariate scatter plots as AMALGAMPar['d'] = 1")
        else:
            print(f"AMALGAM WARNING: Cannot plot map with Pareto correlation estimates as AMALGAMPar['d'] = {AMALGAMPar['d']} [too large]")

        std_X1 = np.std(X1,axis=0)  # Compute the standard deviation of the Pareto solutions
        dx_move = 1/12              # How much do we extend x,y intervals beyond prior range (= 3%)
        # Correlation plots of the posterior parameter samples
        # X,Y axis of scatter plots correspond to prior ranges if defined
        if AMALGAMPar['d'] <= 25:
            # Create a matrix plot for the marginal distributions and bivariate scatter plots
            fig, axs = plt.subplots(AMALGAMPar['d'], AMALGAMPar['d'], figsize=(15, 15))
            # Calculate the number of bins for each parameter
            Nbins = np.array([calcnbins(X1[:, i]) for i in range(AMALGAMPar['d'])])
            nbins = min(np.min(Nbins), maxbins)
            nbins = max(5, round(nbins / 2))  # Adjust number of bins
            for i in range(AMALGAMPar['d']):
                for j in range(AMALGAMPar['d']):
                    if i != j:
                        axs[i, j].scatter(X1[:, j], X1[:, i], color='gray', s=12)
                        # Add a least-squares line for off-diagonal plots
                        # You can use numpy's polyfit for fitting a line if necessary
                        if std_X1[i] > 0 and std_X1[j] > 0:         # only if std exceeds zero
                            fit = np.polyfit(X1[:, j], X1[:, i], 1)
                            axs[i, j].plot(X1[:, j], np.polyval(fit, X1[:, j]), 'r--', linewidth=1)
                        # Add prior ranges if hypercube is specified
                        if 'boundhandling' in Par_info:
                            # Plot vertical lines for parameter bounds
                            axs[i, j].plot([Par_info['min'][j], Par_info['min'][j]], [Par_info['min'][i], Par_info['max'][i]], color='blue', linestyle='-', linewidth=1)
                            axs[i, j].plot([Par_info['max'][j], Par_info['max'][j]], [Par_info['min'][i], Par_info['max'][i]], color='blue', linestyle='-', linewidth=1)
                            axs[i, j].plot([Par_info['min'][j], Par_info['max'][j]], [Par_info['min'][i], Par_info['min'][i]], color='blue', linestyle='-', linewidth=1)
                            axs[i, j].plot([Par_info['min'][j], Par_info['max'][j]], [Par_info['max'][i], Par_info['max'][i]], color='blue', linestyle='-', linewidth=1)
                            # adjust the axes
                            xlow = Par_info['min'] - dx_move * (Par_info['max'] - Par_info['min'])
                            xup = Par_info['max'] + dx_move * (Par_info['max'] - Par_info['min'])
                            axs[i, j].set_xlim([xlow[j], xup[j]])
                            axs[i, j].set_ylim([xlow[i], xup[i]])
                        else:
                            # adjust the axes
                            axs[i, j].set_xlim([min(X1[:, j]), max(X1[:, j])])
                            axs[i, j].set_ylim([min(X1[:, i]), max(X1[:, i])])
                    if i == j:
                        # make a histogram
                        axs[i, j].hist(X1[:, i], bins=nbins, density=True, alpha=0.6, color='gray', edgecolor='black')
                        # Add prior ranges if hypercube is specified
                        if 'boundhandling' in Par_info:
                            # adjust the axes
                            xlow = Par_info['min'] - dx_move * (Par_info['max'] - Par_info['min'])
                            xup = Par_info['max'] + dx_move * (Par_info['max'] - Par_info['min'])
                            axs[i, j].set_xlim([xlow[i], xup[i]])
                        else:
                            # adjust the axes
                            axs[i, j].set_xlim([min(X1[:, i]), max(X1[:, i])])
                    # Set custom x-ticks and x-tick labels
                    x_min, x_max = axs[i, j].get_xlim()
                    dx = x_max - x_min
                    xticks = np.array([x_min + 1/12*dx, x_min + 6/12*dx, x_min + 11/12*dx])
                    axs[i, j].set_xticks(xticks)  # Set x-tick positions first, then labels, otherwise warning
                    axs[i, j].set_xticklabels([str(round(tick,2)) for tick in xticks])
                    y_min, y_max = axs[i, j].get_ylim()
                    dy = y_max - y_min
                    yticks = np.array([y_min + 1/12*dy, y_min + 6/12*dy, y_min + 11/12*dy])
                    axs[i, j].set_yticks(yticks)  # Set y-tick positions first, then labels, otherwise warning
                    axs[i, j].set_yticklabels([str(round(tick,2)) for tick in yticks])
                    # Add values and labels to the axes
                    if i == AMALGAMPar['d'] - 1:
                        axs[i, j].set_xlabel(str_par[j], fontsize=fontsize_xylabel)
                    else:
                        axs[i, j].set_xticklabels([])

                    # Add values and labels to the axes
                    if j == 0:
                        axs[i, j].set_ylabel(str_par[i], fontsize=fontsize_xylabel)
                    else:
                        axs[i, j].set_yticklabels([])

            # Title of the figure
            fig.suptitle(f"AMALGAM: Marginal distribution and bivariate scatter plots of Pareto samples", fontsize=fontsize_title, y=0.92)
            pdf.savefig()
            plt.show()
#            plt.close()
        else:
            print(f"\nAMALGAM WARNING: Cannot plot bivariate scatter plots as AMALGAMPar['d'] = {AMALGAMPar['d']} (= too large)\n")


        ### Plot convergence to Pareto distribution: IGD statistic
        if len(Fpareto) > 0:
            plt.figure(figsize=(15, 6))
            plt.semilogy(output['IGD'][:, 0], output['IGD'][:, 1], 'r')
            plt.xlabel('Number of generations', fontsize=fontsize_xylabel)
            plt.ylabel('Inverted Generational Distance, IGD', fontsize=fontsize_xylabel)
            plt.title('AMALGAM: Evolution of IGD statistic/diagnostic', fontsize=fontsize_title)
            plt.xlim(0, AMALGAMPar['T'])             
            # Set x-axis to show integer ticks
            plt.gca().xaxis.set_major_locator(MaxNLocator(integer=True))
            plt.grid(True)
            plt.tight_layout()
            pdf.savefig()        
            plt.show()

        ### Plot the normalized Pareto parameter solutions
        fig, ax = plt.subplots(figsize=(15, 8))
        Nrank = X1.shape[0]  # Number of rank 1 solutions
        X1_n = (X1 - Par_info['min']) / (Par_info['max'] - Par_info['min'])  # Normalize the parameters
        # Locate the best individual solutions for each objective
        r = np.nan * np.ones(AMALGAMPar['m'])
        for j in range(AMALGAMPar['m']):
            i = np.argmin(FX1[:, j])  # Find the index of the best solution
            r[j] = i
        # Add rank 1 solutions
        for j in range(Nrank):
            if j == 0:
               ax.plot(np.arange(1, AMALGAMPar['d'] + 1), X1_n[j,:], '-', color='gray', linewidth=1, label='Rank 1',zorder=2)
            else:
               ax.plot(np.arange(1, AMALGAMPar['d'] + 1), X1_n[j,:], '-', color='gray', linewidth=1,zorder=1) 
        # Add Edge solutions in color
        symbols = ['r', 'b', 'g', 'c', 'k', 'm', 'y']
        for j in range(AMALGAMPar['m']):
            j_plus_1 = j+1
            ax.plot(np.arange(1, AMALGAMPar['d'] + 1), X1_n[r[j].astype(int), :], symbols[j], linewidth=2, label=r'$F_{{{}, \text{{opt}}}}$'.format(j_plus_1),zorder=3)
        # Add vertical lines and parameter labels
        for j in range(AMALGAMPar['d']):
            ax.axvline(x=j + 1, color='k', linestyle='--', linewidth=1)
            ax.text(j + 1, -0.06, str_par[j], fontsize=fontsize_text, ha='center', va='bottom', rotation=90,zorder=1)
        # Adjust axis limits
        ax.axis([0.8, AMALGAMPar['d'] +.2, 0, 1])
        ax.set_xticks([])
        ax.set_xticklabels([])
        plt.ylabel('Normalized ranges', fontsize=fontsize_xylabel)
        plt.title('AMALGAM: Normalized ranges of Pareto solution samples', fontsize=fontsize_title)
        # ax.legend(loc='upper right', fontsize=fontsize_legend)
        # Match text color to line color
        legend_handles = [Line2D([0], [0], color='gray', lw=4, label=f'Rank 1')]
        legend_handles[1:AMALGAMPar['m']] = [Line2D([0], [0], color=symbols[i], lw=4, label=r'$F_{{{}, \text{{opt}}}}$'.format(i+1)) for i in range(0,AMALGAMPar['m'])]
        legend = ax.legend(handles=legend_handles, fontsize=fontsize_legend, loc='upper right')
        # Match text color to line color
        symbols_ext = ['gray'] + symbols
        for z, text in enumerate(legend.get_texts()):
            text.set_color(symbols_ext[z])  
        plt.tight_layout()
        pdf.savefig()
        plt.show()

        ### Plot the contribution of each algorithm (selection probabilities)
        plt.figure(figsize=(15, 6))
        ha = plt.plot(output['p_rm'][:, 0], output['p_rm'][:, 1:], linewidth=2)
        # Set the colors
        color_order = plt.cm.tab10.colors
        for i, line in enumerate(ha):
            line.set_color(color_order[i % len(color_order)])
        M = np.max(output['p_rm'][:, 1:], axis=0)  # Max selection probability
        plt.axis([0, output['p_rm'][-1, 0], 0, min(1, 1.1 * np.max(M))])
        plt.xlabel('Number of generations', fontsize=fontsize_xylabel)
        plt.ylabel('Selection probability', fontsize=fontsize_xylabel)
        # Set x-axis to show integer ticks
        plt.gca().xaxis.set_major_locator(MaxNLocator(integer=True))
        plt.title('AMALGAM: Selection probability of each recombination method', fontsize=fontsize_title)
        # Create custom legend handles with colored lines
        legend_handles = [Line2D([0], [0], color=color_order[i % len(color_order)], lw=4, label=f"Method: {AMALGAMPar['rec_methods'][i].upper()}") for i in range(output['p_rm'].shape[1] - 1)]
        legend = plt.legend(handles=legend_handles, fontsize=fontsize_legend, loc='best')
        # Match text color to line color
        for z, text in enumerate(legend.get_texts()):
            text.set_color(color_order[z % len(color_order)])  
        
        plt.grid(True)
        plt.tight_layout()
        pdf.savefig()        
        plt.show()

        ### Plot Pareto simulation uncertainty [if model returns a time series]
        if YX is not None:  # Check if YX is not empty
            sim_min = np.min(YX, axis=0)
            sim_max = np.max(YX, axis=0)
            # Open new figure
            plt.figure(figsize=(15, 6))
            plt.fill_between(range(0, YX.shape[1]), sim_min, sim_max, color=(0.75, 0.75, 0.75), label="Pareto simulation uncertainty")
            # Set axis ranges
            delta_y = np.max(sim_max) - np.min(sim_min)
            if np.min(sim_min) > 10:
                plt.axis([0, YX.shape[1], np.min(sim_min) - 0.1 * delta_y, np.max(sim_max) + 0.1 * delta_y])
            else:
                plt.axis([0, YX.shape[1], 0, np.max(sim_max) + 0.1 * delta_y])
            plt.xlabel(r'Time, $t$', fontsize=fontsize_xylabel)
            plt.ylabel(r'Variable of interest, $y$', fontsize=fontsize_xylabel)
            plt.title('AMALGAM: Pareto simulation uncertainty ranges', fontsize=fontsize_title)
            # Adding the legend
            plt.legend(fontsize=fontsize_legend, loc='upper right')
            pdf.savefig()            
            plt.show()

        # Plot the distribution of each parameter
        id = rank_Z(Z, AMALGAMPar, options)
        XZ1 = Z[id, :AMALGAMPar['d']]

        ### Marginal distributions of Pareto parameter (rank 1) solutions [includes past generations]
        row, col = 2, 4
        idx_y_label = np.arange(0, AMALGAMPar['d'] + 1 , col)
        # Calculate the number of bins for each parameter
        Nbins = [calcnbins(XZ1[:, i]) for i in range(AMALGAMPar['d'])]
        nbins = min(np.min(Nbins), maxbins)
        nbins = max(5, round(nbins / 2))  # Adjust number of bins
        counter = 0
        while counter <= AMALGAMPar['d'] - 1:
            # Create a row x col grid of subplots
            fig, axs = plt.subplots(row, col, figsize=(15, 10))
            # Adjust the spacing between subplots
            plt.subplots_adjust(wspace=0.4, hspace=0.4)  # Increase the space horizontally and vertically
            # Loop over each subplot
            for ax in axs.flat:
                # Compute histogram for parameter j
                M, edges = np.histogram(XZ1[:, counter], bins=nbins, density=True)
                X = 0.5 * (edges[:-1] + edges[1:])
                ax.bar(X, M / np.max(M), width=(edges[1] - edges[0]), color='gray', edgecolor='w', alpha=0.7, zorder=1)
                # Add xlabels = parameter name in Latex
                ax.set_xlabel(str_par[counter], fontsize=fontsize_xylabel, labelpad = 10)
                # Adjust axis limits (tight)
                yticks = np.arange(0, 1.02, 0.2)
                ax.set_yticks(yticks)  # Set x-tick positions
                # ax.set_yticklabels([str_par(tick) for tick in yticks])
                ax.set_yticklabels([str(round(tick,1)) for tick in yticks])
                ax.set_ylim(0, 1.02)                # Adjust y limits with padding
                # Add y-label
                if counter in idx_y_label:
                    ax.set_ylabel('Empirical density', fontsize=fontsize_xylabel, labelpad=10)
                # Plot the MAP value
                for ii in range(AMALGAMPar['m']):
                    ax.plot(MAP[ii, counter], 0, f'{symbol[ii]}x', markersize=12, markeredgewidth=3, linewidth=3,zorder=2, clip_on=False)
                # Add letter
                label_plot = get_label(counter)  # This converts 0 -> 'A', 1 -> 'B', etc.
                ax.text(0.02,0.97, f'({label_plot})', transform=ax.transAxes, fontsize=fontsize_text, horizontalalignment='left', va='top', ha='left')
                counter += 1
                if counter == AMALGAMPar['d']:
                    # delete empty subplots
                    for ax in axs.flat:
                        if not ax.has_data():       # If the axis has no data, remove it
                            fig.delaxes(ax)
                    break
            # Set the title of the figure
            fig.suptitle(r"$\text{AMALGAM: Marginal parameter distributions of Pareto (= rank 1) solutions}$", fontsize=fontsize_title, y=0.98)
            # Optionally adjust spacing to avoid overlap with subplots
            fig.tight_layout() #rect=[0, 0, 1, 0.95])
            pdf.savefig()
            plt.show()
            # plt.close()

    # Open the PDF file
    # plt.open(file_name)

    # Open the final PDFs
    os.startfile(file_name)


# Latin Hypercube Sampling function
def LH_sampling(mn, mx, N):
    ## ################################################################################## ##
    ## Latin hypercube sampling of initial chain states DREAM Package                     ##
    ##                                                                                    ##
    ## SYNOPSIS: R = LH_sampling(mn,mx,N)                                                 ##
    ##  where                                                                             ##
    ##   mn        [input] REQUIRED: 1 x d vector of lower bound values                   ##
    ##   mx        [input] REQUIRED: 1 x d vector of upper bound values                   ##
    ##   N         [input] REQUIRED: # of Latin hypercube samples                         ##
    ##   r         [outpt] Nxd matrix of Latin hypercube samples                          ##
    ##                                                                                    ##
    ## Implementation based on the following paper                                        ##
    ##  Minasny, B., and A.B. McBratney (2006), A conditioned Latin hypercube method for  ##
    ##      sampling in the presence of ancillary information, Computers & Geosciences,   ##
    ##      Vol. 32 (9), pp. 1378-138                                                     ##
    ##                                                                                    ##
    ## © Written by Jasper A. Vrugt, Feb 2007                                             ##
    ## Los Alamos National Laboratory                                                     ##
    ##                                                                                    ##
    ## ################################################################################## ##

    #d = len(mn)                                                
    d = mn.shape[1]                                             # Number of parameters
    rng = np.array(mx) - np.array(mn)                           # 1 x d vector with parameter ranges
    y =  np.random.rand(N, d)                                   # N x d matrix with uniform random labels
    # Python: important change below so that X in bound
    # as list is from 0 - N-1 rather than 1 to N
    # id_ = np.argsort(np.random.rand(N, d), axis=0)            # Sorting random values to avoid duplicates
    id_matrix = 1 + np.argsort(np.random.rand(N, d), axis=0)    # Random sort (1:N without replacement)
    M = (id_matrix - y) / N                                     # Multiplier matrix (y introduces randomness)
    R = np.add(np.multiply(M, rng), mn)  # N x d matrix of stratified LH samples

    return R


def get_label(counter):
    # ####################################################################### #
    # This function turns an integer into a letter/letters                    #
    # ####################################################################### #

    label = ""
    while counter >= 0:
        label = chr(65 + (counter % 26)) + label
        counter = counter // 26 - 1
    return label


def convert_memoryview_to_array(plugin):
    # ####################################################################### #
    # Convert memoryview objects in plugin to numpy arrays to make picklable  #
    # ####################################################################### #

    for key, value in plugin.items():
        if isinstance(value, memoryview):
            plugin[key] = np.array(value)  # Convert memoryview to ndarray
    
    return plugin


def distance(FQ, RQ, options):
    # ####################################################################### #
    # Computes distance (crowding distance or strength Pareto) for points     #
    # ####################################################################### #

    M, m = FQ.shape
    if options['density'] == 'crowding':
        dQ = np.full(M, np.nan)
        for j in range(1, np.max(RQ) + 1):
            id_rank_j = (RQ == j)
            FQ_rank_j = FQ[id_rank_j]
            N_s = np.sum(id_rank_j)
            Cr_d = np.zeros(N_s)
            for i in range(m):
                sorted_values = np.sort(FQ_rank_j[:, i])
                sorted_idx = np.argsort(FQ_rank_j[:, i])
                Cr_d[sorted_idx[[0, N_s - 1]]] = np.inf
                for z in range(1, N_s - 1):
                    Cr_d[sorted_idx[z]] += sorted_values[z + 1] - sorted_values[z - 1]
            dQ[id_rank_j] = Cr_d
    elif options['density'] == 'strength':
        n_dom = np.zeros(M)
        for j in range(M):
            n_dom[j] = np.sum(np.all(FQ[j] < FQ, axis=1))
        dQ = M / n_dom  # Higher n_dom means lower strength (more crowded)
    else:
        raise ValueError("Distance: Unknown density estimation method")

    return dQ


def Update_PS(AMALGAMPar, Z, PS, FZ_min):
    # ####################################################################### #
    # Updates dictionary of Particle Swarm                                    #
    # ####################################################################### #

    # Euclidean distance to FZ_min for each point in Z (considering objective values starting from AMALGAMPar['d'])
    Ed = np.sum((Z[:, AMALGAMPar['d']:AMALGAMPar['d'] + AMALGAMPar['m']] - FZ_min) ** 2, axis=1)

    # Find the minimum Euclidean distance
    sorted_ids = np.argsort(Ed)
    
    # Now take the AMALGAMPar['N'] best values of Ed as individual best each parent
    PS['n'] = Z[sorted_ids[:AMALGAMPar['N']], :AMALGAMPar['d']]

    # Now take the minimum of Ed as overall best solution
    PS['p'] = np.tile(Z[sorted_ids[0], :AMALGAMPar['d']], (AMALGAMPar['N'], 1))

    return PS


def checkfile_T_AMALGAM(fname):
    # ####################################################################### #
    # Check the content of the restart budget file 'T.txt'                    #
    # ####################################################################### #

    # Check if file exists
    if not os.path.exists(fname):
        raise FileNotFoundError(f"AMALGAM ERROR: File '{fname}' does not exist.")

    # Read the content of the file
    with open(fname, 'r') as f:
        content = f.read().strip()

    # Check if the content is empty
    if not content:
        raise ValueError(f"AMALGAM ERROR: File '{fname}' is empty --> Please store an integer in file '{fname}'.")

    # Try to convert the content to a numeric value
    try:
        T_new = int(content)
    except ValueError:
        raise ValueError(f"AMALGAM ERROR: File '{fname}' does not store a numerical value --> Store only a single value (integer) in file '{fname}'.")

    # Check if the value is a positive integer
    if T_new <= 0:
        raise ValueError(f"AMALGAM ERROR: File '{fname}' stores negative integers (or zero) --> Store a single positive integer in file '{fname}'.")

    return T_new


def NSGA(AMALGAMPar, Par_info, X, FX, id, dX):

    G = NSGA_crossover(AMALGAMPar, Par_info, X, FX, dX)     # Selection & crossover
    G = NSGA_mutate(AMALGAMPar, Par_info, G)                # Polynomial mutation
    G = G[id, :AMALGAMPar['d']]                             # Select individuals
    
    return G


def NSGA_crossover(AMALGAMPar, Par_info, X, FX, dX):
    # Perform genetic selection and crossover
    
    div = AMALGAMPar['N'] // 2                              # Divide population size by 2
    res = np.ceil(div - np.floor(div))                      # Handle uneven population size
    
    a1 = np.random.permutation(AMALGAMPar['N'])             # Random permutation for parents
    a2 = np.random.permutation(AMALGAMPar['N'])
    
    # Adjust indices if population size is odd
    a1 = np.concatenate([a1, a1[:int(res)]])
    a2 = np.concatenate([a2, a2[:int(res)]])
    
    G = []                                                  # List to store offspring
    ct = 0
    
    # Perform tournament selection and crossover
    while len(G) < AMALGAMPar['N']:
        for j in range(2):
            if j == 0:
                a_1, a_2, a_3, a_4 = a1[ct], a1[ct + 1], a1[ct + 2], a1[ct + 3]
            else:
                a_1, a_2, a_3, a_4 = a2[ct], a2[ct + 1], a2[ct + 2], a2[ct + 3]
            
            parent1 = Tournament(AMALGAMPar, X[a_1], X[a_2], FX[a_1], FX[a_2], dX[a_1], dX[a_2])
            parent2 = Tournament(AMALGAMPar, X[a_3], X[a_4], FX[a_3], FX[a_4], dX[a_3], dX[a_4])
            
            # Perform crossover between parents
            child1, child2 = crossover(AMALGAMPar, Par_info, parent1, parent2)
            
            # Add children to offspring
            G.extend([child1, child2])
            ct += 1
    
    # Select first N individuals
    G = np.array(G)[:AMALGAMPar['N'], :AMALGAMPar['d']]

    return G


def Tournament(AMALGAMPar, a, b, Fa, Fb, da, db):
    # Perform tournament selection
    
    rnd = np.random.rand()
    done = False
    
    flagout = Check_dominance(AMALGAMPar, Fa, Fb)
    
    if flagout == 1:
        parent = a
        done = True
    elif flagout == -1:
        parent = b
        done = True
    
    if not done:
        if da > db:
            parent = a
            done = True
        elif da < db:
            parent = b
            done = True
    
    if not done:
        if rnd <= 0.5:
            parent = a
        else:
            parent = b
    
    return parent


def Check_dominance(AMALGAMPar, Fa, Fb):
    # Check for dominance in multi-objective sense
    
    flag1, flag2 = 0, 0
    
    for qq in range(AMALGAMPar['m']):
        if Fa[qq] < Fb[qq]:
            flag1 = 1
        elif Fa[qq] > Fb[qq]:
            flag2 = 1
    
    delta = flag1 - flag2
    flagout = -1
    
    if delta == 0:
        flagout = 0
    elif delta == 1:
        flagout = 1
    
    return flagout


def crossover(AMALGAMPar, Par_info, x1, x2):
    # Perform crossover between two individuals
    
    if np.random.rand() < AMALGAMPar['p_CR']:
        yl = Par_info['min'][0, :AMALGAMPar['d']]
        yu = Par_info['max'][0, :AMALGAMPar['d']]
        
        rnd = np.random.rand(AMALGAMPar['d'])
        idx = np.where(rnd > 0.5)[0]
        idx_par = np.where(np.abs(x1 - x2) <= np.finfo(float).eps)[0]
        
        y1 = np.minimum(x1, x2)
        y2 = np.maximum(x1, x2)
 
         # adjusted by JAV to avoid dividing by zero
        # beta = 1.0 + (2.0 * (y1 - yl) / (y2 - y1))
        beta = 1.0 + (2.0 * (y1 - yl) / ((y2 - y1) + 1e-10) )

        alpha = 2.0 - np.power(beta, -(AMALGAMPar['eta_C'] + 1.0))
        
        rnd = np.random.rand(AMALGAMPar['d'])
        ii1 = np.where(rnd <= (1.0 / alpha))[0]
        ii2 = np.where(rnd > (1.0 / alpha))[0]
        
        betaq = np.ones(AMALGAMPar['d']) * 1e-10
        betaq[ii1] = np.power(rnd[ii1] * alpha[ii1], 1.0 / (AMALGAMPar['eta_C'] + 1.0))
        betaq[ii2] = np.power(1.0 / (2.0 - rnd[ii2] * alpha[ii2]), 1.0 / (AMALGAMPar['eta_C'] + 1.0))
        
        g1 = 0.5 * ((y1 + y2) - betaq * (y2 - y1))
        
        # adjusted by JAV to avoid dividing by zero
#        beta = 1.0 + (2.0 * (yu - y2) / (y2 - y1))
        beta = 1.0 + (2.0 * (yu - y2) / ((y2 - y1) + 1e-10) )
        alpha = 2.0 - np.power(beta, -(AMALGAMPar['eta_C'] + 1.0))
        
        betaq[ii1] = np.power(rnd[ii1] * alpha[ii1], 1.0 / (AMALGAMPar['eta_C'] + 1.0))
        betaq[ii2] = np.power(1.0 / (2.0 - rnd[ii2] * alpha[ii2]), 1.0 / (AMALGAMPar['eta_C'] + 1.0))
        
        g2 = 0.5 * ((y1 + y2) + betaq * (y2 - y1))
        
        # Bound children
        g1 = np.maximum(g1, yl)
        g1 = np.minimum(g1, yu)
        g2 = np.maximum(g2, yl)
        g2 = np.minimum(g2, yu)
        
        rnd = np.random.rand(AMALGAMPar['d'])
        
        ii1 = np.where(rnd <= 0.5)[0]
        ii2 = np.where(rnd > 0.5)[0]
        
        g1[ii1] = g2[ii1]
        g2[ii1] = g1[ii1]
        
#        g1[ii2] = g1[ii2]
#        g2[ii2] = g2[ii2]
        
        g1[idx_par] = x1[idx_par]
        g2[idx_par] = x2[idx_par]
        
        g1[idx] = x1[idx]
        g2[idx] = x2[idx]

    else:
        g1 = x1
        g2 = x2
    
    return g1, g2


def NSGA_mutate(AMALGAMPar, Par_info, G):
    # Perform polynomial mutation
    
    delta1 = (G - Par_info['min']) / (Par_info['max'] - Par_info['min'])
    delta2 = (Par_info['max'] - G) / (Par_info['max'] - Par_info['min'])
    
    rnd = np.random.rand(AMALGAMPar['N'], AMALGAMPar['d'])
    mut_pow = 1 / (AMALGAMPar['eta_M'] + 1)

    # Added by JAV
    xy = np.zeros(G.shape)
    val = np.zeros(G.shape)
    deltaq = np.zeros(G.shape)

    idx = np.where(rnd <= 0.5)
#    xy = 1 - delta1[idx]
    xy[idx] = 1 - delta1[idx]
#    val = 2 * rnd[idx] + (1 - 2 * rnd[idx]) * (xy ** (AMALGAMPar['eta_M'] + 1))
    val[idx] = 2 * rnd[idx] + (1 - 2 * rnd[idx]) * (xy[idx] ** (AMALGAMPar['eta_M'] + 1))
#    deltaq = val ** mut_pow - 1
    deltaq[idx] = val[idx] ** mut_pow - 1
    
    idx = np.where(rnd > 0.5)
#    xy = 1 - delta2[idx]
    xy[idx] = 1 - delta2[idx]
#    val = 2 * (1 - rnd[idx]) + 2 * (rnd[idx] - 0.5) * (xy ** (AMALGAMPar['eta_M'] + 1))
    val[idx] = 2 * (1 - rnd[idx]) + 2 * (rnd[idx] - 0.5) * (xy[idx] ** (AMALGAMPar['eta_M'] + 1))
#    deltaq = 1 - val ** mut_pow
    deltaq[idx] = 1 - val[idx] ** mut_pow
    
    z = G + deltaq.reshape(AMALGAMPar['N'], AMALGAMPar['d']) * (Par_info['max'] - Par_info['min'])
    
    rnd = np.random.rand(AMALGAMPar['N'], AMALGAMPar['d'])
    idx = np.where(rnd <= AMALGAMPar['p_M'])
    
    G[idx] = z[idx]
    
    return G


def AM(AMALGAMPar, X, R, n, id):
    """
    Generate offspring using the adaptive Metropolis algorithm.
    
    Parameters:
    AMALGAMPar : object or dict with the necessary parameters
    X : numpy.ndarray
        Current population (matrix where rows are individuals and columns are variables).
    R : numpy.ndarray
        A selection mask, typically indicating which individuals to use for covariance estimation.
    n : int
        The number of offspring to generate.
    id : numpy.ndarray
        Indices of the selected individuals for mutation.
    
    Returns:
    G : numpy.ndarray
        The generated offspring.
    """
    # Calculate the covariance structure using selected individuals
    selected_X = X[R == 1, :AMALGAMPar['d']]            # Select individuals for covariance estimation
    if selected_X.shape[0] > 1:                         # Make sure there is more than 1 row
        cov_matrix = np.cov(selected_X.T)
        cov_matrix += 1e-10 * np.eye(AMALGAMPar['d'])   # Add a small value to avoid singular matrix
        R_matrix = np.linalg.cholesky(cov_matrix)       # Perform Cholesky decomposition
    else:
        R_matrix = np.eye(AMALGAMPar['d'])
    
    # Generate children
    G = X[id, :AMALGAMPar['d']] + np.random.randn(n, AMALGAMPar['d']) @ (AMALGAMPar['gamma'] * R_matrix)

    return G


def PSO(AMALGAMPar, X, PS, n, id):
    """
    Generate offspring using Particle Swarm Optimization (PSO) method.
    
    Parameters:
    AMALGAMPar : object or dict with the necessary parameters
    X : numpy.ndarray
        Current population (matrix where rows are individuals and columns are variables).
    PS : object or dict with particle swarm states
        Contains velocity (v), best positions (p and n).
    n : int
        The number of offspring to generate.
    id : numpy.ndarray
        Indices of the selected individuals for mutation.
    
    Returns:
    G : numpy.ndarray
        The generated offspring.
    PS : object or dict
        The updated particle swarm states.
    """
    # Draw random numbers
    r_1 = np.random.rand(AMALGAMPar['N'], AMALGAMPar['d'])
    r_2 = np.random.rand(AMALGAMPar['N'], AMALGAMPar['d'])
    # Draw varphi (scaling factors)
    # varphi = AMALGAMPar['varphi'](AMALGAMPar['N'])
    varphi = np.random.uniform(AMALGAMPar['varphi']-1/4, AMALGAMPar['varphi']+1/4, AMALGAMPar['N']).reshape(-1,1) 
    
    # Compute velocity of the swarm
    PS['v'] = np.multiply(PS['v'], varphi) + AMALGAMPar['c_1'] * r_1 * (PS['p'] - X) + AMALGAMPar['c_2'] * r_2 * (PS['n'] - X)
    # Update particle positions
    G = X[id, :AMALGAMPar['d']] + PS['v'][id, :AMALGAMPar['d']]
    
    # Apply random rotation
    xi = 1 + np.random.uniform(-1, 1, (n, 1))  # Random values between -1 and 1
    G = np.multiply(G, xi)  # Apply element-wise multiplication

    return G, PS


def DE(AMALGAMPar, X, n, id):
    """
    Generate offspring using Differential Evolution (DE) method.

    Parameters:
    AMALGAMPar : object or dict with the necessary parameters
        Contains `beta_1`, `beta_2`, and `d` values.
    X : numpy.ndarray
        Current population (matrix where rows are individuals and columns are variables).
    n : int
        The number of offspring to generate.
    id : numpy.ndarray
        Indices of the selected individuals for mutation.

    Returns:
    G : numpy.ndarray
        The generated offspring.
    """
    #F = AMALGAMPar['beta_1'](n) * np.ones((1, AMALGAMPar['d']))
    F = np.random.uniform(AMALGAMPar['beta_1'] - 0.2,AMALGAMPar['beta_1'] + 0.2, n).reshape(-1,1) * np.ones((1, AMALGAMPar['d'])),
    #K = AMALGAMPar['beta_2'](n) * np.ones((1, AMALGAMPar['d']))
    K = np.random.uniform(AMALGAMPar['beta_2'] - 0.2,AMALGAMPar['beta_2'] + 0.2, n).reshape(-1,1) * np.ones((1, AMALGAMPar['d'])),
    # Generate random values and sort to create rr
    rr = np.argsort(np.random.rand(n, AMALGAMPar['N']), axis=1)
     # Generate children using differential evolution
    G = X[id, :AMALGAMPar['d']] + F * (X[rr[:, 0], :AMALGAMPar['d']] - X[id, :AMALGAMPar['d']]) + \
        K * (X[rr[:, 1], :AMALGAMPar['d']] - X[rr[:, 2], :AMALGAMPar['d']])

    return G


def Boundary_handling(X, Par_info):
    # ####################################################################### #
    # This function checks that parameters are in prior bounds, corrects them #
    # ####################################################################### #
    """
    Parameters:
    X (numpy.ndarray): Required, N x d matrix of candidate points
    Par_info (dict): Required, dictionary containing parameter bounds (min/max) and boundary treatment method
    
    Returns:
    Xr (numpy.ndarray): N x d matrix with revised candidate points
    v (numpy.ndarray): N x 1 vector with 0 for in bound and 1 for out of bound
    """
    
    Xr = np.copy(X)                                         # Create a copy of X for the revised values
    N, d = X.shape
    v = np.zeros(N, dtype=bool)                             # Logical array indicating out-of-bound values
    
    mn = np.tile(Par_info['min'], (N, 1))                   # Lower bounds replicated N times
    mx = np.tile(Par_info['max'], (N, 1))                   # Upper bounds replicated N times
    
    # Positions where X is below lower bound or above upper bound
    id_l = np.where(X < mn)                                 # Smaller than lower bound
    id_u = np.where(X > mx)                                 # Larger than upper bound
    
    # Boundary handling options
    if Par_info['boundhandling'] == 'reflect':              # Reflection method
        Xr[id_l] = 2 * mn[id_l] - X[id_l]                   # Reflect below the lower bound
        Xr[id_u] = 2 * mx[id_u] - X[id_u]                   # Reflect above the upper bound
    elif Par_info['boundhandling'] == 'bound':              # Bound method
        Xr[id_l] = mn[id_l]                                 # Set to lower bound
        Xr[id_u] = mx[id_u]                                 # Set to upper bound
    elif Par_info['boundhandling'] == 'fold':               # Folding method
        Xr[id_l] = mx[id_l] - (mn[id_l] - X[id_l])          # Fold below the lower bound
        Xr[id_u] = mn[id_u] + (X[id_u] - mx[id_u])          # Fold above the upper bound
    elif Par_info['boundhandling'] == 'reject':             # Reject method
        o = np.zeros_like(X, dtype=bool)                    # Initialize out-of-bound array
        o[id_l] = 1                                         # Mark positions below the lower bound
        o[id_u] = 1                                         # Mark positions above the upper bound
        v = np.sum(o, axis = 1) > 0                           # Identify rows with any out-of-bound values
    
    # Reflection or folding: Check if all elements are within bounds
    # Both methods can go out of bounds if violation exceeds |mx - mn|
    if Par_info['boundhandling'] in ['reflect', 'fold']:
        id_l = np.where(Xr < mn)                            # Smaller than lower bound
        id_u = np.where(Xr > mx)                            # Larger than upper bound
        Xr[id_l] = np.random.uniform(mn[id_l], mx[id_l])    # Random draw in [mn, mx]
        Xr[id_u] = np.random.uniform(mn[id_u], mx[id_u])    # Random draw in [mn, mx]

    # BMA model training if applicable
    # if 'unit_simplex' in Par_info:
    #     wght_sum = np.sum(Xr[:N, :int(K)], axis = 1)
    #     Xr[:, :int(K)] = Xr[:, :int(K)] / wght_sum[:, np.newaxis]  # Normalize weights in the unit simplex

    return Xr, v


def plot_screen(fig, ax0, AMALGAMPar, F, t, F_par, output):
    """
    Plot AMALGAM results to screen during each trial.
    
    Parameters:
    fig : matplotlib.figure.Figure
        The figure object.
    ax0 : matplotlib.axes.Axes
        The axes for the plot.
    AMALGAMPar : dict
        Contains algorithm parameters such as 'rec_methods' and 'q'.
    F : numpy.ndarray
        The current objective values.
    t : int
        The current generation/trial.
    F_par : numpy.ndarray
        Pareto front values for comparison.
    output : dict
        Contains algorithm outputs like 'IGD' and 'p_alg'.
        
    Returns:
    fig, ax0 : matplotlib.figure.Figure, matplotlib.axes.Axes
        Updated figure and axes.
    """
    
    names = ['GA', 'PSO', 'AMS', 'DE']              # Individual methods' names
    symbol = ['rs', 'gd', 'mp', 'co', 'bh']         # Marker types for different methods
    Colors = np.array([[1, 0, 0], [0, 0.5, 0], [0.5, 0.5, 0.5], [1, 0.64, 0], [0, 0, 1]])
    
    gif = True  # Enable GIF output
    
    # Identify the index of the selected method
    if len(AMALGAMPar['rec_methods']) == 1:
        method = AMALGAMPar['rec_methods'][0]
        if method == 'ga':
            idx = 0
        elif method == 'pso':
            idx = 1
        elif method == 'ams':
            idx = 2
        elif method == 'de':
            idx = 3
    else:
        idx = 4  # AMALGAM
        method = 'AMALGAM'
    
    if t == 2:
        fig, ax0 = plot_ZDT(F_par)
        if gif:
            capture_frame(fig, method)
    
    ax0 = plt.gca()  # Use the current axes

    # Plot the results for the current method
    ax0.semilogy(F[:, 0], F[:, 1], symbol[idx], color = Colors[idx], linewidth = 1, markersize = 8)
    ax0.set_xlim(0, 1)
    ax0.set_ylim(1e-4, 1e3)
    
    IGD_print = f"{round(10000 * output['IGD'][t, 2]) / 10000:.4f}"
    IGD_text = ax0.text(0.85, 1600, f"${IGD_print}$", fontsize=16, ha='center')
    
    ax0.text(0.65, 1600, f"${t}$", fontsize=16, ha='center')

    # Progress bar
    ax1 = fig.add_axes([0.05, 0.7, 0.1, 0.05])  # Position of the bar
    ax1.barh(0, t / AMALGAMPar['T'], color=[0.25, 0.25, 0.25])
    ax1.set_xlim(0, 1)
    ax1.set_yticks([])
    ax1.set_xticks([])

    # Probability bar for individual methods
    p_alg = output['p_alg'][t, 1:AMALGAMPar['q']+1]
    if idx < 4:
        P_alg = np.zeros(4)
        P_alg[idx] = p_alg
    else:
        P_alg = p_alg
    
    # Plot individual selection probabilities for methods
    for qq in range(4):
        ax_pos = [0.1, 0.55 - qq * 0.1, 0.8, 0.05]
        ax1 = fig.add_axes(ax_pos)
        ax1.barh(0, P_alg[qq], color=Colors[qq])
        ax1.set_xlim(0, 1)
        ax1.set_yticks([])
        ax1.set_xticks([])
        ax1.set_xticklabels([0, 0.2, 0.4, 0.6, 0.8, 1.0])
        if qq == 3:
            ax1.set_xticklabels([f'{i:.1f}' for i in np.linspace(0, 1, 6)])
            ax1.text(-0.2, 0, 'SELECTION PROBABILITY', fontsize=14)

        ax1.text(-0.2, 0, names[qq], color=Colors[qq], rotation=0, fontweight='bold', fontsize=14)

    if gif:
        capture_frame(fig, method, append=True)

    plt.draw()

    # Clear previous plots for the next iteration
    IGD_text.remove()
    
    return fig, ax0


def capture_frame(fig, method, append=False):
    """
    Capture the current frame and save as a GIF.
    """
    frame = fig.canvas.draw()
    img = np.array(fig.canvas.renderer.buffer_rgba())
    img = img[:, :, :3]
    img_ind, cmap = plt.matplotlib.colors.rgb_to_hsv(img)
    
    mode = 'append' if append else 'w'
    plt.imsave(f'{method}_ZDT.gif', img_ind, cmap=cmap, format='gif', mode=mode)


def plot_ZDT(F_par):
    """
    Plot the true Pareto front for ZDT function.
    """
    fig, ax = plt.subplots(figsize=(6, 9))
    ax.semilogy(F_par[:, 0], F_par[:, 1], 'ko', markersize=2, linewidth=2, markeredgecolor='k', markerfacecolor='k')
    ax.set_xlim(0, 1)
    ax.set_ylim(1e-4, 1e3)
    ax.set_xlabel('$F_1(\theta)$', fontsize=18)
    ax.set_ylabel('$F_2(\theta)$', fontsize=18)
    ax.tick_params(axis='both', which='major', labelsize=16)
    
    return fig, ax


def prt_progress(AMALGAMPar, N, count=[0], ct=[0]):
    """
    Print the progress of model simulations.
    
    Parameters:
    AMALGAMPar : dict
        Contains the parameters, specifically the CPU usage for progress tracking.
    N : int
        Total number of simulations.
    count : list
        Persistent counter for printing the progress. Default is [0].
    ct : list
        Persistent counter for the current count of simulations. Default is [0].
    """
    # Update the count and ct (simulating persistent behavior)
    ct[0] += AMALGAMPar['CPU']
    
    # Handle potential issue where ct doesn't reach N because of remaining balance
    # (This is commented in the original MATLAB code and isn't implemented here, as it's a more complex case)
    
    # Clear the line before printing the new progress (simulates the \b behavior)
    print('\r' + ' ' * count[0], end='')  # '\r' moves the cursor to the beginning of the line
    
    # Update count for how many characters have been printed (for backspace handling)
    count[0] = print(f'Model simulations, %% done: {100 * (ct[0] / N):3.2f}', end='', flush=True)


def Discrete_space(X, Par_info, method = 2):
    """
    This function transforms continuous space of X to discrete values.

    Parameters:
    X : numpy.ndarray
        A N x d matrix of candidate points (continuous space).
    Par_info : dict
        A dictionary containing parameter information, including:
        - 'min': The lower bounds for each dimension.
        - 'max': The upper bounds for each dimension.
        - 'steps': The number of steps for discretization.
        - 'step_size': The size of each step in the discretized space.
    
    Returns:
    X_dis : numpy.ndarray
        A N x d matrix of discrete parameter values.
    """
    
    N = X.shape[0]  # Number of candidate vectors

    # Step 1: Transform continuous X to integer between 0 and number of steps
    # Step 2: Back transform to discrete space
    
    if method == 1:
        # Proper method for all MATLAB releases (older version)
        X_min = np.tile(Par_info['min'], (N, 1))  # Replicate min values N times
        X_int = np.round(np.tile(Par_info['steps'], (N, 1)) * ((X - X_min) / np.tile(Par_info['max'] - Par_info['min'], (N, 1))))
        X_dis = X_min + X_int * np.tile(Par_info['step_size'], (N, 1))
    
    elif method == 2:
        # New MATLAB method (for later releases)
        X_int = np.round(Par_info['steps'] * ((X - Par_info['min']) / (Par_info['max'] - Par_info['min'])))
        X_dis = Par_info['min'] + X_int * Par_info['step_size']
    
    return X_dis


def calcnbins(x, method='middle', minb=1, maxb=np.inf):
    """
    Compute the "ideal" number of bins using different methods for histogram bin calculation.
    
    Parameters:
    - x: vector of data points
    - method: string with choice of method for calculating bins. Default is 'middle'.
        Options: 'fd', 'scott', 'sturges', 'middle', 'all'
    - minb: smallest acceptable number of bins (default: 1)
    - maxb: largest acceptable number of bins (default: np.inf)
    
    Returns:
    - nbins: The calculated number of bins based on the selected method.
    """
    
    # Input checking
    if not isinstance(x, (np.ndarray, list, np.generic)):
        raise ValueError('The x argument must be numeric or logical.')

    x = np.asarray(x)
    
    # Ensure the array is real, discard imaginary part
    if np.iscomplexobj(x):
        x = np.real(x)
        print('Warning: Imaginary parts of x will be ignored.')
    
    # Ensure x is a vector (1D array)
    if x.ndim != 1:
        x = x.flatten()
        print('Warning: x will be coerced to a vector.')
    
    # Remove NaN values
    x = x[~np.isnan(x)]
    if len(x) == 0:
        raise ValueError("x must contain at least one valid number.")
    
    # Choose method if not specified
    valid_methods = ['fd', 'scott', 'sturges', 'all', 'middle']
    if method not in valid_methods:
        raise ValueError(f"Unknown method: {method}")
    
    # Method selection
    if method == 'fd':
        nbins = calc_fd(x)
    elif method == 'scott':
        nbins = calc_scott(x)
    elif method == 'sturges':
        nbins = calc_sturges(x)
    elif method == 'middle':
        nbins = [calc_fd(x), calc_scott(x), calc_sturges(x)]
        nbins = np.median(nbins)
    elif method == 'all':
        nbins = {'fd': calc_fd(x),
                'scott': calc_scott(x),
                'sturges': calc_sturges(x)}
    
    # Confine number of bins to the acceptable range
    nbins = confine_to_range(nbins, minb, maxb)
    
    return nbins


def calc_fd(x):
    """Freedman-Diaconis rule"""
    h = np.subtract(*np.percentile(x, [75, 25]))  # Interquartile range (IQR)
    if h == 0:
        h = 2 * np.median(np.abs(x - np.median(x)))  # Median absolute deviation (MAD)
    
    if h > 0:
        nbins = np.ceil((np.max(x) - np.min(x)) / (2 * h * len(x) ** (-1/3)))
    else:
        nbins = 1
    return nbins


def calc_scott(x):
    """Scott's method"""
    h = 3.5 * np.std(x) * len(x) ** (-1/3)
    if h > 0:
        nbins = np.ceil((np.max(x) - np.min(x)) / h)
    else:
        nbins = 1
    return nbins


def calc_sturges(x):
    """Sturges' method"""
    nbins = np.ceil(np.log2(len(x)) + 1)
    return nbins


def confine_to_range(x, lower, upper):
    """Ensure bin count is within the specified range"""
    x = np.maximum(x, lower)
    x = np.minimum(x, upper)
    return np.floor(x)


def safe_int(value):
    return int(value) if value else 0  # or return a default value if empty


def safe_float(value):
    return float(value) if value else 0.0  # or return a default value if empty


def rank_Z(Z, AMALGAMPar, options):
    """
    rank_Z: Function which ranks all past populations iteratively
    Written by Jasper A. Vrugt
    University of California Irvine
    """

    N, d2 = Z.shape  # Determine size of Z
    d = d2 - AMALGAMPar['m']  # # parameters
    m = AMALGAMPar['m']  # # objective functions
    Nrank = 1000  # Maximum # points to rank at once

    if N <= Nrank:
        # Rank Z at once
        R = AMALGAM_rank(Z, options)
        id = np.where(R == 1)[0]
    else:
        # Take segments of Z
        q = N // Nrank
        R1 = np.zeros(N, dtype=bool)  # Create a boolean array for ranks
        for ii in range(q):
            id_s = ii * Nrank
            id_e = (ii + 1) * Nrank
            R, _, _ = AMALGAM_rank(Z[id_s:id_e, d:d + m], options)
            R1[id_s:id_e] = (R == 1)

        if id_e < N:
            R, _, _ = AMALGAM_rank(Z[id_e:N, d:d + m], options)
            R1[id_e:N] = (R == 1)

        # Locate indices of rank 1 solutions
        id_1 = np.where(R1 == 1)[0]

        # Now rank all rank 1 solutions of different segments of Z at once
        R, _, _ = AMALGAM_rank(Z[id_1, d:d + m], options)
        
        # Which are now the rank 1 solutions of Z?
        id_1F = (R == 1)
        id = id_1[id_1F]

    return id


def IGD(S, Q, C = None):
    """
    Calculates the Inverse Generational Distance (IGD) metric for a set of solutions.

    Args:
        S (numpy.ndarray): The ideal Pareto front (PF*), shape (row, Scol)
        Q (numpy.ndarray): The obtained nondominated front (PF), shape (row, Qcol)
        C (numpy.ndarray, optional): Constraint matrix, shape (Crow, Qcol). Default is None.

    Returns:
        float: The IGD value
    """
    row, Scol = S.shape  # S is the ideal Pareto front (PF*)
    row, Qcol = Q.shape  # Q is the obtained front (PF)

    # Step 1: Remove infeasible solutions based on the constraints
    feasible = np.ones(Qcol, dtype=int)     # All solutions are initially feasible
    nof = Qcol                              # Number of feasible solutions
    if C is not None:
        for i in range(Qcol):
            for k in range(C.shape[0]):     # Iterate over rows of C
                if C[k, i] < -1E-6:         # Infeasible solution
                    feasible[i] = 0         # Mark as infeasible
                    nof -= 1                # Decrease the count of feasible solutions
                    break
    if nof == 0:
        return 1.0E6  # Return large value if no feasible solutions

    # Step 2: Calculate the IGD value for feasible solutions
    dis = 0.0
    for i in range(Scol):                   # For each solution in the ideal Pareto front
        min_dist = 1.0E200                  # Start with a large number
        for j in range(Qcol):               # Compare to each solution in the obtained front
            if feasible[j] > 0:             # Only consider feasible solutions
                # Calculate the Euclidean distance
                d = np.sum((S[:, i] - Q[:, j])**2)
                min_dist = min(min_dist, d)  # Keep the minimum distance
        dis += np.sqrt(min_dist)  # Add the minimum distance to the total

    # Normalize by the number of solutions in S
    return dis / Scol


def Compute_FX_true(AMALGAMPar, Func_name, problem, M=500, plugin=None):
    """
    Computes actual Pareto optimal solutions.

    Args:
        AMALGAMPar (dict): Structure with AMALGAM settings/parameters.
        Func_name (str): Name of the function that computes objective functions.
        problem (str): Name of benchmark problem, such as 'ZDT1', 'ZDT2', etc.
        M (int, optional): Number of points on the true Pareto front (default is 500).
        plugin (object, optional): 2nd argument of Func_name: Class set by user.

    Returns:
        FX_true (ndarray): True Pareto front (approximation) in objective space.
    """

    if plugin is None:
        plugin = {}

    options = {'restart': 'no',
                'parallel': 'no',
                'modout': 'no'}
    
    M = int(M)
    AMALGAMPar['d'] = int(AMALGAMPar['d'])
    AMALGAMPar, func_handle, base_dir = AMALGAM_calc_setup(AMALGAMPar, Func_name, options, plugin)
    AMALGAMPar['N'] = int(AMALGAMPar['N'])
    AMALGAMPar['m'] = int(AMALGAMPar['m']) 
    AMALGAMPar['T'] = int(AMALGAMPar['T'])
    N_max = (10000)                             # Number of points to sample
    
    # Define minimum values for FX based on the problem
    if problem.upper() in ['ZDT1', 'ZDT2', 'ZDT3', 'ZDT4', 'ZDT6']:     ## Zitzler's problems
        #FX1_min = np.zeros(1, AMALGAMPar['d'])
        FX1_min = np.zeros(AMALGAMPar['d']).reshape(1, AMALGAMPar['d'])
        FX2_min = np.array([1] + [0] * (AMALGAMPar['d'] - 1)).reshape(1, AMALGAMPar['d'])
    elif problem.upper() in ['DTLZ1', 'DTLZ2', 'DTLZ3']:                ## Deb's problems
        FX1_min = np.concatenate(([0, 0], 0.5 * np.ones(AMALGAMPar['d'] - 2))).reshape(1, AMALGAMPar['d'])
        FX2_min = np.concatenate(([1, 1], 0.5 * np.ones(AMALGAMPar['d'] - 2))).reshape(1, AMALGAMPar['d'])
    else:
        raise ValueError("Unknown problem")

    T = FX1_min.shape[0]
    N = [0] + [N_max // T] * T                          # Number of samples
    Ncs = np.cumsum(N)                                  # Cumulative sum of samples
    X = np.full((Ncs[-1], AMALGAMPar['d']), np.nan)     # Initialize parameter values
    
    for j in range(0,T):
        # Latin hypercube sampling
        Xlh = LH_sampling(FX1_min[j,:].reshape(1, -1), FX2_min[j].reshape(1, -1), N[j+1])
        Xlh[0, :AMALGAMPar['d']] = FX1_min[j, :AMALGAMPar['d']]    # Single objective optimum
        Xlh[1, :AMALGAMPar['d']] = FX2_min[j, :AMALGAMPar['d']]    # Single objective optimum
        X[Ncs[j]:Ncs[j+1], :AMALGAMPar['d']] = Xlh                 # Store samples
 
    printed_warnings = set()
    # Compute objective values (you need to define AMALGAM_calc_FX)
    FX = AMALGAM_calc_FX(X, AMALGAMPar, options, func_handle, base_dir, plugin, printed_warnings)[0]

    # Get true Pareto front approximation (you need to define getF_true)
    FX_true = getF_true(FX, AMALGAMPar, 1, M)
    
    return FX_true


def getF_true(FX, AMALGAMPar, N_f, M):
    """
    Generates uniformly spaced Pareto solutions.
    
    Args:
        FX (ndarray): Sorted objective space values.
        AMALGAMPar (dict): AMALGAM parameters.
        N_f (int): Number of fronts to approximate.
        M (int): Number of points for Pareto approximation.
    
    Returns:
        FX_true (ndarray): True Pareto front (approximated).
    """

    m = AMALGAMPar['m']             # Number of objectives
    FX = FX[FX[:, 0].argsort()]     # Sort FX by the first objective

    # Find dissimilar FX points
    id = np.where(np.diff(FX[:, 0]) > 0)[0]
    FX = np.vstack([FX[0, :m], FX[id + 1, :m]])

    N_FX = FX.shape[0]              # Number of dissimilar points
    Ed = np.zeros(N_FX)             # Euclidean distance between adjacent points
    Ed[0] = 0
    for i in range(1, N_FX):
        Ed[i] = np.linalg.norm(FX[i, :m] - FX[i - 1, :m])

    # Sort distances
    Ed_s, ii_s = np.sort(-Ed), np.argsort(-Ed)

    # Determine front to sample
    if N_f == 1:
        id = [0, N_FX]
        Ed_cs = np.cumsum(Ed)
    else:
        id = [0] + list(ii_s[:N_f - 1]) + [N_FX + 1]
        Ed_cs = np.cumsum(-Ed_s[N_f:N_FX])

    mEd = Ed_cs[-1] / (M - 1)           # Mean distance between points
    FX_true = np.full((M, m), np.nan)   # Initialize Pareto approximation   
    nq = np.zeros(N_f + 1).astype(int)
    for i in range(N_f):
        Ed_f = np.cumsum(Ed[id[i] + 1:id[i + 1] - 1])
        Ed_cs = Ed_f[-1]
        # Define FX of the front
        FX_f = FX[id[i] + 1:id[i + 1] - 1, :m]
        nf = FX_f.shape[0]  # Number of points
        if mEd > Ed_cs:
            F_pareto = FX_f[0, :]  # For interpolation
        else:
            Ed_int = np.arange(mEd, Ed_cs, mEd)
            #F_pareto = np.full((Ed_int.shape[0], m), np.nan)
            F_pareto = np.full((M-2, m), np.nan)
            for z in range(m):
                F_pareto[:,z] = np.interp(Ed_int, Ed_f, FX_f[:,z])

        F_pareto = np.vstack([FX_f[0, :m], F_pareto])   # First point F_pareto
        F_pareto = np.vstack([F_pareto, FX_f[-1, :m]])  # Last point F_pareto
        q = ~np.isnan(F_pareto[:, 0])
        nq[i + 1] = np.sum(q)
        FX_true[nq[i]:nq[i + 1], :] = F_pareto[q,:m]
    
    return FX_true


def distribute_tasks(N, CPU):
    # ####################################################################### #
    # Split the task ranges evenly across workers                             #
    # ####################################################################### #

    chunk_size = N // CPU
    task_ranges = []
    for i in range(CPU):
        start_idx = i * chunk_size
        # Ensure the last worker gets the remaining tasks
        end_idx = N if i == CPU - 1 else (i + 1) * chunk_size
        task_ranges.append((start_idx, end_idx))
    return task_ranges


############################### Interface with user-defined function ############################### 

def get_function_handle(Func_name):
    """
    This function returns a function handle (reference to the function)
    for the specified Func_name. This avoids repeated imports and lookups.
    """
    # Check if Func_name contains a dot (.) to separate module and function   
    if '.' in Func_name:
        module_name, function_name = Func_name.rsplit('.', 1)
    else:
        # If no dot is present, assume Func_name is the function name in the current module
        module_name = __name__      # Use the current module
        function_name = Func_name
    
    try:
        # Dynamically import the module if it's not already imported
        module = importlib.import_module(module_name)
    except ModuleNotFoundError:
        raise ValueError(f"Module '{module_name}' not found.")
    
    # Get the function from the module
    Func = getattr(module, function_name, None)
    
    if Func is None:
        raise ValueError(f"Function '{function_name}' not found in module '{module_name}'.")
    
    # Ensure the function is callable
    if not callable(Func):
        raise ValueError(f"The object '{function_name}' in '{module_name}' is not callable.")
    
    # Return the function handle (reference to the function)
    return Func


def worker(ii, func_handle, X, plugin=None):

    if plugin is not None:
        result = func_handle(X[ii, :], plugin)
    else:
        result = func_handle(X[ii, :])

    return result

# Function to recursively convert memoryview objects to bytearrays
def convert_memoryview(obj):
    if isinstance(obj, memoryview):
        # Convert memoryview to bytearray
        return bytearray(obj)
    elif isinstance(obj, dict):
        # Recursively convert memoryview objects inside a dictionary
        return {key: convert_memoryview(value) for key, value in obj.items()}
    elif isinstance(obj, list):
        # Recursively convert memoryview objects inside a list
        return [convert_memoryview(item) for item in obj]
    elif isinstance(obj, tuple):
        # Recursively convert memoryview objects inside a tuple
        return tuple(convert_memoryview(item) for item in obj)
    elif isinstance(obj, np.ndarray):
        # Convert memoryview from NumPy arrays to ndarray (if necessary)
        return np.array(obj)
    else:
        # Return other objects as they are (assumed to be picklable)
        return obj


def copy_model_files(source_dir, target_dir):
    """
    Copy model files from the source directory to the target worker directory.
    """
    if not os.path.exists(target_dir):
        os.makedirs(target_dir)
    
    # List all files in the source directory and copy them to the target directory
    for filename in os.listdir(source_dir):
        source_file = os.path.join(source_dir, filename)
        target_file = os.path.join(target_dir, filename)
        
        if os.path.isfile(source_file):
            shutil.copy2(source_file, target_file)  # Copy file with metadata
        elif os.path.isdir(source_file):
            shutil.copytree(source_file, target_file)  # Copy directory recursively


def worker_task(worker_id, start_idx, end_idx, X, func_handle, plugin, base_dir):
    
    if base_dir is not None:
        # Create a worker-specific directory
        worker_dir = os.path.join(base_dir, f"worker_{worker_id}")   
        # Change to the worker's specific directory
        os.chdir(worker_dir)

    # Execute model
    results = []
    for idx in range(start_idx, end_idx):
        if plugin is not None:
            result = func_handle(X[idx, :], plugin)
        else:
            result = func_handle(X[idx, :])

        results.append(result)

    # Return the results for this worker
    return results


def cleanup_worker_directories(base_dir, N):
    """
    Clean up (delete) the worker directories after all generations are complete.
    """
    for worker_id in range(N):
        worker_dir = os.path.join(base_dir, f'worker_{worker_id}')
        # Delete the worker directory
        if os.path.exists(worker_dir):
            shutil.rmtree(worker_dir)  


############################# End Interface with user-defined function ############################# 
