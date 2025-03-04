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
#  SYNOPSIS                                                               #
#   [X,FX,output,Z,YX] = AMALGAM(AMALGAMPar,Func_name,Par_info);          #
#   [X,FX,output,Z,YX] = AMALGAM(AMALGAMPar,Func_name,Par_info, ...       #
#       options);                                                         #
#   [X,FX,output,Z,YX] = AMALGAM(AMALGAMPar,Func_name,Par_info, ...       #
#       options,plugin);                                                  #
#   [X,FX,output,Z,YX] = AMALGAM(AMALGAMPar,Func_name,Par_info, ...       #
#       options,plugin,Ftrue);                                            #
#  where                                                                  #
#   AMALGAMPar  [input] Structure with AMALGAM settings/parameters        #
#    .d             Number of parameters                                  #
#    .T             Number of generations                                 #
#    .m             Number of objective functions                         #
#    .N             Population size                                       #
#    .rec_methods   Recombination methods                                 #
#      = 'ga'       Genetic algorithm                                     #
#      = 'am'       Adaptive Metropolis                                   #
#      = 'de'       Differential evolution                                #
#      = 'ps'       Particle swarm optimization                           #
#    .beta1         Scaling factor for DE algorithm       DEF: U[0.6,1]   #
#    .beta2         Scaling factor for DE algorithm       DEF: U[0.2,0.6] #
#    .c1            Social factor, PSO algorithm          DEF: 1.5        #
#    .c2            Cognitive factor, PSO algorithm       DEF: 1.5        #
#    .varphi        Inertia factor, PSO algorithm         DEF: U[0.5,1.0] #
#    .p_CR          Crossover probability, GA algorithm   DEF: 0.9        #
#    .p_M           Mutation probability, GA algorithm    DEF: 1/d        #
#    .eta_C         Crossover index, GA algorithm         DEF: 10         #
#    .eta_M         Mutation index, GA algorithm          DEF: 50         #
#    .gamma         Jumprate AMS algorithm                DEF: 2.38^2/d   #
#    .K             Thinning constant                     DEF: 1          #
#    .p0            Minimum probability crossover methods DEF: 0.05       #
#   Func_name   [input] Name (string) function computes objective functs. #
#   Par_info    [input] Parameter structure: Ranges, initial/prior & bnd  #
#    .names         1xd-cell array with parameter names   DEF: []         #
#    .min           1xd-vector of min parameter values    DEF: -inf(1,d)  #
#    .max           1xd-vector of max parameter values    DEF: inf(1,d)   #
#    .boundhandling Treat the parameter bounds or not?                    #
#      = 'reflect'  Reflection method                                     #
#      = 'bound'    Set to bound                                          #
#      = 'fold'     Folding [Vrugt&Braak: doi:10.5194/hess-15-3701-2011]  #
#      = 'none'     No boundary handling                  DEFault         #
#    .initial       Method to draw initial chain states                   #
#      = 'uniform'  Uniform: U(Par_info.min,Par_info.max)                 #
#      = 'latin'    Latin hypercube: LH(Par_info.min,Par_info.max)        #
#      = 'normal'   Normal:  N(Par_info.mu,Par_info.cov)                  #
#      = 'prior'    User specified prior distribution                     #
#      = 'user'     Initial population taken from Par_info.x0             #
#    .mu            1xd-mean vector: µ if .initial = 'normal'             #
#    .cov           dxd-covariance matrix: Σ if .initial = 'normal'       #
#    .prior         Prior distribution (manual) if .initial = 'prior'     #
#                   Ex 1: Par_info.prior = @(x,a,b) mvnpdf(x,a,b);        #
#                             Par_info.a = [-2 -2];                       #
#                             Par_info.b = eye(2);                        #
#                   Ex 2: Par_info.prior = {'normpdf(x,0,1)',...          #
#                                           'unifpdf(x,-2,2)'}            #
#    .steps         d-vector with # intervals for each parameter          #
#                    → Discrete AMALGAM                                   #
#   options     [input] Structure with computational settings/options     #
#    .parallel      Multi-core computation chains?        DEF: 'yes'      #
#    .IO            If parallel, IO writing model?        DEF: 'no'       #
#    .screen        Print screen output during trial?     DEF: 'no'       #
#    .density       Density of rank 1,2,3, etc. solutions DEF: 'crowding' #
#     = 'crowding'  Crowding distance:Deb et al: NSGA-II  DEFault         #
#     = 'strength'  Strength Pareto:Zitzler&Thiele: SPEA-2                #
#    .modout        Return model simulations?             DEF: 'no'       #
#    .save          Save AMALGAM output during the run?   DEF: 'no'       #
#    .restart       Restart run? (only with "save")       DEF: 'no'       #
#    .print         Output writing screen (tables/figs)   DEF: 'yes'      #
#   plugin      [input] OPT: 2nd argument of Func_name: Class set by user #
#   Ftrue       [input] OPT: NxM matrix true Pareto front, objectve space #
#   X           [outpt] Nxd matrix final population                       #
#   FX          [outpt] Nxm matrix final objective functions              #
#   output      [outpt] Structure summarizes algorithmic performance      #
#    .p_alg         Selection probability crossover methods               #
#    .IGD           Inverse generational distance [if F_par defined]      #
#    .RunTime   [outpt] CPU time in seconds                               #
#   Z           [outpt] NxTxd matrix of past populations                  #
#   YX          [outpt] Model simulations of Pareto solutions             #
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
import os, sys
import time

example_dir = os.getcwd()					                    # Add current directory to Python path
if example_dir not in sys.path:
    sys.path.append(example_dir)

parent_dir = os.path.abspath(os.path.join(example_dir, '..'))   # Go up one directory
sys.path.append(os.path.join(parent_dir, 'miscellaneous'))	    # Add miscellaneous directory to Python path
from AMALGAM_functions import *				                    # Import functions


#def AMALGAM(AMALGAMPar, Func_name, Par_info, *args):
    # Initialize options and plugin as empty
#    options = {'restart': 'no'}
#    plugin = {}
#    Ftrue = []
#    file_name = 'AMALGAM'

    # Handle input arguments
#    if len(args) == 1:
#        options = args[0]
#    elif len(args) == 2:
#        options, plugin = args
#    elif len(args) == 3:
#        options, plugin, Ftrue = args

#    if 'restart' not in options:
#        options['restart'] = 'no'

## Main program
def AMALGAM(AMALGAMPar, Func_name, Par_info, options = None, plugin = None, Ftrue = []):
    # ####################################################################### #
    # Main program (function) of AMALGAM - Python implementation              #
    # ####################################################################### #

    if options == [] or options == None:
        options = {}
    if 'restart' not in options:
        options['restart'] = 'no'
    if plugin == []:
        plugin = None
    if Ftrue is None:
        Ftrue = []    

    printed_warnings = set()                                                                                # Initialize a set to track warnings for Evaluate_target
    file_name = 'AMALGAM'                                                                                   # Name of memory dump
    if options['restart'] == 'no':                                                                          # New trial: No restart

        AMALGAMPar, Par_info, options = AMALGAM_check(Func_name, AMALGAMPar, Par_info, options, Ftrue)      # Check input variables
        AMALGAMPar, Par_info, options, T_start = AMALGAM_setup(AMALGAMPar, Par_info, options)               # Define all algorithmic variables
        AMALGAMPar, func_handle, base_dir = AMALGAM_calc_setup(AMALGAMPar, Func_name, options, plugin)      # Initialize computational environment    
        AMALGAMPar, Par_info, X, p_rm, PS, Z, output = AMALGAM_initialize(AMALGAMPar, Par_info, plugin)     # Initialize all variables
        FX, YX = AMALGAM_calc_FX(X, AMALGAMPar, options, func_handle, base_dir, plugin, printed_warnings)   # Compute objective functions initial population
        RX, dX, FX_min = AMALGAM_rank(FX, options)                                                          # Rank initial population
        Z[:AMALGAMPar['N'], :AMALGAMPar['d'] + AMALGAMPar['m']] = np.concatenate([X, FX], axis=1)           # Store initial population in archive
        
        if 'ps' in AMALGAMPar['rec_methods']:                                                               
            PS = Update_PS(AMALGAMPar, Z[:AMALGAMPar['N'], :AMALGAMPar['d'] + AMALGAMPar['m']], PS, FX_min) # Update dictionary of Particle Swarm

        if len(Ftrue) > 0:  
            output['IGD'][0, :] = np.concatenate([np.array([0]), np.array([IGD(Ftrue.T, FX.T)])])           # Compute hypervolume

        ct = 1                                                                                              # Counter for external archive of generations    
    elif options['restart'] == 'yes':                                                                       # Restart run [= continue where stopped]
        AMALGAMPar, Par_info, func_handle, options, PS, X, Z, FX, YX, p_rm, output, FX_min, RX, \
            dX, ct, base_dir, T_start = AMALGAM_restart(file_name, Func_name, Ftrue)

    t0 = time.time()
    fg, ax0 = None, None  # for plotting purposes
 
    # Dynamic Part
    for t in range(T_start, AMALGAMPar['T'] + 1):                                                           # Iterating through generations
        id, id_rm = AMALGAM_distribution(AMALGAMPar, p_rm)                                                  # Create distribution rec. methods
        G, PS = AMALGAM_children(AMALGAMPar, Par_info, X, FX, RX, dX, PS, id_rm, plugin)                    # Create children with different recombination methods
        FG, YG = AMALGAM_calc_FX(G, AMALGAMPar, options, func_handle, base_dir, plugin, printed_warnings)   # Compute objective functions of children
        FXG_min = np.minimum(FX_min, np.min(FG, axis=0))                                                    # Minimum values of each objective function [= minimization]
        X, FX, RX, dX, id_N, id = AMALGAM_population(AMALGAMPar, options, X, G, FX, FG, id)                 # New population 
        p_rm = AMALGAM_load(AMALGAMPar, p_rm, id)                                                           # New selection probability rec. methods
        
        if options.get('modout', 'no') == 'yes':
            simtot = np.vstack([YX, YG])
            YX = simtot[id_N, :]

        if t % AMALGAMPar['K'] == 0:                                                                        # Append to archive
            Z[ct * AMALGAMPar['N']: (ct + 1) * AMALGAMPar['N'], :AMALGAMPar['d'] + AMALGAMPar['m']] \
                = np.concatenate([X, FX], axis=1)
            ct += 1
            if 'ps' in AMALGAMPar['rec_methods']:                                                           # Update Particle Swarm dictionary
                PS = Update_PS(AMALGAMPar, Z[:ct * AMALGAMPar['N'], :AMALGAMPar['d'] + AMALGAMPar['m']], PS, FXG_min)
        output['p_rm'][t, 0:AMALGAMPar['q'] + 1] = np.concatenate([[t], p_rm]) 
        
        if len(Ftrue) > 0:                                                                                  # Compute hypervolume
            output['IGD'][t, :] = np.concatenate([np.array([t]), np.array([IGD(Ftrue.T, FX.T)])])  
        
        if options['save'] == 'yes':                # Open a shelve file to store the data                  # Save output to file
            np.save(file_name, {                    # with shelve.open(file_name, 'c') as file:
                    'AMALGAMPar': AMALGAMPar,               # file['AMALGAMPar'] = AMALGAMPar
                    'Par_info': Par_info,                   # file['Par_info'] = Par_info
                    'options': options,                     # file['options'] = options
                    'PS': PS,                               # file['PS'] = PS
                    'output': output,                       # file['output'] = output
                    'X': X,                                 # file['X'] = X
                    'Z': Z,                                 # file['Z'] = Z
                    'FX': FX,                               # file['FX'] = FX
                    'YX': YX,                               # file['YX'] = YX
                    'FX_min': FX_min,                       # file['FX_min'] = FX_min
                    'RX': RX,                               # file['RX'] = RX
                    'dX': dX,                               # file['dX'] = dX
                    'p_rm': p_rm,                           # file['p_rm'] = p_rm
                    't': t,                                 # file['t'] = t
                    'ct': ct})                              # file['ct'] = ct

        if options.get('screen', 'no') == 'yes':                                                            # Plot to screen [= animation]
            fg, ax0 = plot_screen(fg, ax0, AMALGAMPar, FX, t, Ftrue, output)
        
        if t % (max(1,AMALGAMPar["T"] // 25) ) == 0:                                                        # Print progress
            if t > 0:
                print(f'AMALGAM calculating, {100*(t/AMALGAMPar["T"]):.2f}% done', end='\r')

    print('\n')
    output['RunTime'] = time.time() - t0                                                                    # Determine total run time
    AMALGAM_end(AMALGAMPar, options, base_dir)                                                              # Close workers
    Z = Z[:ct * AMALGAMPar['N'], :AMALGAMPar['d'] + AMALGAMPar['m']]                                        # Finalize external archive

    if options['print'] == 'yes':                                                                           # Print progress of postprocessing
        for t in range(2, 4):
            if t == 2:
                print_name = '........'
            if t == 3:
                print_name = '........ done'
            if t > 1:
                print(f"\r\nAMALGAM postprocessor, please wait {print_name}", end="\n")
                if t == 2:
                    AMALGAM_postproc(AMALGAMPar, Par_info, options, X, FX, output, Ftrue, YX, Z)            # Create Tables and Figures               

    return X, FX, output, Z, YX
