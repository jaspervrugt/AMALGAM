% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% %
%                                                                         %
%   EEEEEE  XX  XX   AAAA   MM   MM  PPPPPP  LL      EEEEEE   11  33333   %
%   EE       XXXX   AA  AA  MMM MMM  PP  PP  LL      EE       11     33   %
%   EEEEE     XX    AA  AA  MMMMMMM  PPPPPP  LL      EEEEE    11    333   %
%   EE       XXXX   AAAAAA  MM   MM  PP      LL      EE       11     33   %
%   EEEEEE  XX  XX  AA  AA  MM   MM  PP      LLLLLL  EEEEEE   11  33333   %
%                                                                         %
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% %

% Multiple criteria BMA model training: river discharge ensemble
%  Vrugt, J.A. (2024), Distribution-Based Model Evaluation and            
%      Diagnostics: Elicitability, Propriety, and Scoring Rules for       
%      Hydrograph Functionals, Water Resources Research, 60,              
%      e2023WR036710, https://doi.org/10.1029/2023WR036710
%  Vrugt, J.A., and B.A. Robinson (2007), Treatment of uncertainty using 
%      ensemble methods: Comparison of sequential data assimilation and 
%      Bayesian model averaging, Water Resources Research, 43, W01411, 
%      https://doi.org/10.1029/2005WR004838
%  Vrugt, J.A., M.P. Clark, C.G.H. Diks, Q. Duan, and B. A. Robinson      
%      (2006), Multi-objective calibration of forecast ensembles using     
%      Bayesian model averaging, Geophysical Research Letters, 33,        
%      L19817, https://doi.org/10.1029/2006GL027126

clc; clear; close all hidden            % clear memory and figures

AMALGAMPar.N = 100;                     % Define population size
AMALGAMPar.T = 100;                     % # generations?
AMALGAMPar.m = 3;                       % # objective functions?

Par_info.initial = 'latin';             % Latin hypercube sampling
Par_info.boundhandling = 'reflect';     % Explicit boundary handling

Func_name = 'AMALGAM_BMA';              % Define name of function

% Data multi-criteria BMA model training 
data = load('discharge.txt');           % Daily discharge forecasts (mm/d) 
                                        % verifying data 
T_idx = 1:1:3000;                       % Start/end day of training period
D = data(T_idx,1:8);                    % Ensemble forecasts
y = data(T_idx,9);                      % Verifying observations

options.BMA = 'yes';                    % Activate BMA method
PDF = 'normal';                         % Forecast pdf ('normal'/'gamma')
VAR = '4';                              % Variance option ('1'/'2'/'3'/'4')
[AMALGAMPar,Par_info,D_bc,A,B] = ...    % Setup BMA model (+ bias crrction)
    setup_BMA(AMALGAMPar, ...
    Par_info,D,y,VAR);
plugin.BMA = struct('PDF',PDF, ...      % Structure with BMA info AMALGAM
    'VAR',VAR,'D',D_bc,'y',y, ...
    'K',size(D,2));

options.print = 'yes';                  % Print output to screen (figures)
options.modout = 'yes';                 % Return hmodel simulations
options.ranking = 'C';                  % Use Pareto ranking in C
options.save = 'yes';			        % Save memory AMALGAM restart run 	

% Run the AMALGAM code and obtain non-dominated solution set
[X,F,output,Z] = AMALGAM(AMALGAMPar,Func_name,Par_info,options, ...
    plugin);
