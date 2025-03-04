% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% %
%                                                                         %
%   EEEEEE  XX  XX   AAAA   MM   MM  PPPPPP  LL      EEEEEE        4444   %
%   EE       XXXX   AA  AA  MMM MMM  PP  PP  LL      EE           44 44   %
%   EEEEE     XX    AA  AA  MMMMMMM  PPPPPP  LL      EEEEE       44  44   %
%   EE       XXXX   AAAAAA  MM   MM  PP      LL      EE          444444   %
%   EEEEEE  XX  XX  AA  AA  MM   MM  PP      LLLLLL  EEEEEE          44   %
%                                                                         %
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% %

% ZDT3 example from the following paper
%  Zitzler, E., K. Deb, and L. Thiele (2000), Comparison of Multiobjective
%      Evolutionary Algorithms: Empirical Results, Evolutionary 
%      Computation, 8 (2), 183-195, 2000
%      https://sop.tik.ee.ethz.ch/publicationListFiles/zdt2000a.pdf

clc; clear; close all hidden            % clear memory and figures

AMALGAMPar.N = 100;                     % Define population size
AMALGAMPar.T = 100;                     % # generations?
AMALGAMPar.d = 10;                      % # parameters?
AMALGAMPar.m = 2;                       % # objective functions?

Par_info.initial = 'latin';             % Latin hypercube sampling
Par_info.boundhandling = 'bound';       % Explicit boundary handling
Par_info.min = [0 -5*ones(1,AMALGAMPar.d-1)];   % Minimum parameter values
Par_info.max = [1  5*ones(1,AMALGAMPar.d-1)];   % Maximum parameter values

Func_name = 'AMALGAM_ZDT4';             % Define name of function
plugin = AMALGAMPar.d;                  % Plugin: 2nd input arg. Func_name
FX_true = Compute_FX_true( ...          % Create true Pareto front
    AMALGAMPar,Func_name, ...           % = benchmark problem
    'ZDT4',500,plugin);

options.print = 'yes';                  % Print output to screen (figures)
options.ranking = 'C';                  % Use C script for ranking (faster)

% Run the AMALGAM code and obtain non-dominated solution set
[X,FX,output,Z] = AMALGAM(AMALGAMPar,Func_name,Par_info,options, ...
    plugin,FX_true);
