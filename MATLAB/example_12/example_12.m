% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% %
%                                                                         %
%   EEEEEE  XX  XX   AAAA   MM   MM  PPPPPP  LL      EEEEEE   11  222222  %
%   EE       XXXX   AA  AA  MMM MMM  PP  PP  LL      EE       11  22 22   %
%   EEEEE     XX    AA  AA  MMMMMMM  PPPPPP  LL      EEEEE    11    22    %
%   EE       XXXX   AAAAAA  MM   MM  PP      LL      EE       11   22     %
%   EEEEEE  XX  XX  AA  AA  MM   MM  PP      LLLLLL  EEEEEE   11  222222  %
%                                                                         %
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% %

% Multiple criteria hymod model training: driven & nondriven hydrograph
%  Vrugt, J.A., H.V. Gupta, L.A. Bastidas, W. Bouten, and S. Sorooshian   
%      (2000), Effective and efficient algorithm for multiobjective       
%      optimization of hydrologic models, Water Resources Research,       
%      39 (8), 1214, https://doi.org./10.1029/2002WR001746                

clc; clear; close all hidden;           % clear workspace and figures

AMALGAMPar.N = 100;                     % Define population size
AMALGAMPar.T = 150;                     % # generations?
AMALGAMPar.d = 5;                       % # parameters?
AMALGAMPar.m = 2;                       % # objective functions?

Par_info.initial = 'latin';             % Latin hypercube sampling
Par_info.boundhandling = 'reflect';     % Explicit boundary handling
Par_info.min = [1.0 0.10 0.10 0.00 0.10];   % Minimum parameter values
Par_info.max = [500 2.00 0.99 0.10 0.99];   % Maximum parameter values
Par_info.steps = 499*ones(1,AMALGAMPar.d);

Func_name = 'AMALGAM_hymod';            % Define name of function
                                        % Fixed Euler integration
                                        % Better implementations available            
% Data multi-criteria hymod training 
load bound.txt;                         % Daily data Leaf River watershed (7/28/1952 - 9/30/1962)
Area = 1944;                            % Area of Leaf River in km^2
conv_mult = ...                         % Convert discharge m3/s to mm/day 
    Area*(1000*1000)/(1000*60*60*24);
T_max = 795;                            % Use two years of data
Y_obs = bound(65:T_max,4)/conv_mult;    % Measured discarge in mm/day
PET = bound(1:T_max,5);                 % Pot. evptrnspirtion, PET (mm/day)
R = sum(bound(1:T_max,6:9),2);          % Daily rainfall (mm/day) (= Σ 6-hourly data)
idx_d = find(R(65:end) > 0);            % Indexes driven part hydrograph
N_d = numel(idx_d);                     % # entries driven part
idx_nd = find(R(65:end) == 0);          % Indexes nondriven part hydrograph
N_nd = numel(idx_nd);                   % # entries nondriven part
field_names = {'fieldnames', ...        % Collect fields of plugin
    'T_max','Y_obs','PET','R', ...
    'idx_d','N_d','idx_nd','N_nd'};
plugin = v2struct(T_max,Y_obs,PET, ...  % Define plugin
    R,idx_d,N_d,idx_nd,N_nd,field_names); 
plugin.fields = field_names;            % fieldnames of plugin

options.print = 'yes';                  % Print output to screen (figures)
options.parallel = 'yes';               % Multi-core: No IO writing
options.IO = 'no';                      % No I/O writing workers (default)
options.modout = 'yes';                 % Return hmodel simulations
options.ranking = 'C';                  % Example: Use Pareto ranking in C

% Run the AMALGAM code and obtain non-dominated solution set
[X,F,output,Z] = AMALGAM(AMALGAMPar,Func_name,Par_info,options, ...
    plugin);
