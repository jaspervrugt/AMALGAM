% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% %
%                                                                         %
%   EEEEEE  XX  XX   AAAA   MM   MM  PPPPPP  LL      EEEEEE      888888   %
%   EE       XXXX   AA  AA  MMM MMM  PP  PP  LL      EE          88  88   %
%   EEEEE     XX    AA  AA  MMMMMMM  PPPPPP  LL      EEEEE       888888   %
%   EE       XXXX   AAAAAA  MM   MM  PP      LL      EE              88   %
%   EEEEEE  XX  XX  AA  AA  MM   MM  PP      LLLLLL  EEEEEE      888888   %
%                                                                         %
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% %

% DTLZ3 example from the following paper
%  Deb, K., L. Thiele, M. Laumanns, and E. Zitzler (2001), Scalable Test  
%      Problems for Evolutionary Multi-Objective Optimization. Kanpur,    
%      India: Kanpur Genetic Algorithms Lab. (KanGAL), Indian Institute   
%      of Technology, KanGAL Report 2001001                                

clc; clear; close all hidden            % clear memory and figures

AMALGAMPar.N = 250;                     % Define population size
AMALGAMPar.T = 500;                     % # generations?
AMALGAMPar.d = 10;                      % # parameters?
AMALGAMPar.m = 3;                       % # objective functions?

Par_info.initial = 'latin';             % Latin hypercube sampling
Par_info.boundhandling = 'bound';       % Explicit boundary handling
Par_info.min = zeros(1,AMALGAMPar.d);   % Minimum parameter values
Par_info.max = ones(1,AMALGAMPar.d);    % Maximum parameter values

Func_name = 'AMALGAM_DTLZ3';            % Define name of function
plugin = AMALGAMPar.m;                  % Plugin: 2nd input arg. Func_name
FX_true = Compute_FX_true( ...          % Create true Pareto front
    AMALGAMPar,Func_name, ...           % = benchmark problem
    'DTLZ3',500,plugin);                % 3D problems: not well gridded

options.print = 'yes';                  % Print output to screen (figures)
options.ranking = 'C';                  % Use C script for ranking (faster)

% Run the AMALGAM code and obtain non-dominated solution set
[X,FX,output,Z] = AMALGAM(AMALGAMPar,Func_name,Par_info,options, ...
    plugin,FX_true);
