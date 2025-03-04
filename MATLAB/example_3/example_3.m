% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% %
%                                                                         %
%   EEEEEE  XX  XX   AAAA   MM   MM  PPPPPP  LL      EEEEEE      333333   %
%   EE       XXXX   AA  AA  MMM MMM  PP  PP  LL      EE              33   %
%   EEEEE     XX    AA  AA  MMMMMMM  PPPPPP  LL      EEEEE          333   %
%   EE       XXXX   AAAAAA  MM   MM  PP      LL      EE              33   %
%   EEEEEE  XX  XX  AA  AA  MM   MM  PP      LLLLLL  EEEEEE      333333   %
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
AMALGAMPar.d = 30;                      % # parameters?
AMALGAMPar.m = 2;                       % # objective functions?

Par_info.initial = 'latin';             % Latin hypercube sampling
Par_info.boundhandling = 'bound';       % Explicit boundary handling
Par_info.min = zeros(1,AMALGAMPar.d);   % If 'latin', min values
Par_info.max = ones(1,AMALGAMPar.d);    % If 'latin', max values

Func_name = 'AMALGAM_ZDT3';             % Define name of function
plugin = AMALGAMPar.d;                  % Plugin: 2nd input arg. Func_name
FX_true = Compute_FX_true( ...          % Create true Pareto front
    AMALGAMPar,Func_name, ...           % = benchmark problem
    'ZDT3',500,plugin);

options.print = 'yes';                  % Print output to screen (figures)

% Run the AMALGAM code and obtain non-dominated solution set
[X,FX,output,Z] = AMALGAM(AMALGAMPar,Func_name,Par_info,options, ...
    plugin,FX_true);
