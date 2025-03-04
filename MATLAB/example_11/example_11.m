% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% %
%                                                                         %
%   EEEEEE  XX  XX   AAAA   MM   MM  PPPPPP  LL      EEEEEE     11  11    %
%   EE       XXXX   AA  AA  MMM MMM  PP  PP  LL      EE         11  11    %
%   EEEEE     XX    AA  AA  MMMMMMM  PPPPPP  LL      EEEEE      11  11    %
%   EE       XXXX   AAAAAA  MM   MM  PP      LL      EE         11  11    %
%   EEEEEE  XX  XX  AA  AA  MM   MM  PP      LLLLLL  EEEEEE     11  11    %
%                                                                         %
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% %

% Hydrologic modeling using hmodel: Driven & nondriven hydrograph
%  Vrugt, J.A., H.V. Gupta, L.A. Bastidas, W. Bouten, and S. Sorooshian
%      (2000), Effective and efficient algorithm for multiobjective 
%      optimization of hydrologic models, Water Resources Research, 
%      39 (8), 1214, https://doi.org./10.1029/2002WR001746

clc; clear; close all hidden            % clear memory and figures

AMALGAMPar.N = 100;                     % Define population size
AMALGAMPar.T = 150;                     % # generations?
AMALGAMPar.d = 7;                       % # parameters?
AMALGAMPar.m = 2;                       % # objective functions?

Par_info.initial = 'latin';             % Latin hypercube sampling
Par_info.boundhandling = 'reflect';     % Explicit boundary handling
Par_info.min = [1   10  0.1 0.1 -10 0.1 0.1];   % Minimum parameter values
Par_info.max = [10 1000 100 100  10 10  150];   % Maximum parameter values
Par_info.names = {'I_\text{max}', ...
    'S_\text{u,max}','Q_\text{s,max}', ...
    '\alpha_\text{E}','\alpha_\text{F}', ...
    'K_\text{f}','K_\text{s}'};
Func_name = 'AMALGAM_hmodel';           % Define name of function

% Data multi-criteria hmodel training 
mopex = load('03451500.dly');           % Now load the mopex data
idx = find(mopex(:,1) > 1959 & ...
    mopex(:,1) < 1999);                 % Define training data set
n = size(mopex,1); tout = 0:n;          % # observations and time
Y_obs = mopex(idx(1:n),6);              % Measured discharge data (mm/day)
Y_obs = Y_obs(731:n)';                   
data.P     = mopex(idx(1:n),4)';        % Daily rainfall (mm/d)
data.Ep    = mopex(idx(1:n),5)';        % Daily evaporation (mm/d)
data.aS    = 1e-6;                      % Percolation coefficient
hmodel_opt.InitialStep = 1;             % Initial time-step (d)
hmodel_opt.MaxStep     = 1;             % Maximum time-step (d)
hmodel_opt.MinStep     = 1e-6;          % Minimum time-step (d)
hmodel_opt.RelTol      = 1e-3;          % Relative tolerance
hmodel_opt.AbsTol      = 1e-3; %*ones(5,1);% Absolute tolerances (mm)
hmodel_opt.Order       = 2;             % 2nd order method (Heun)
y0 = 1e-6*ones(5,1);                    % Initial conditions
id_d = find(data.P(731:n)>0);           % Index driven part
N_d = numel(id_d);                      % # observations driven part    
id_nd = find(data.P(731:n)==0);         % Index nondriven part
N_nd = numel(id_nd);                    % # observations nondriven part
field_names = {'fieldNames','tout', ... % Create fields of plugin 
    'data','hmodel_opt','y0', ...
    'Y_obs','n','id_d','N_d', ...
    'id_nd','N_nd'};                           
plugin = v2struct(tout,data, ...
    hmodel_opt,y0,Y_obs,n,id_d, ...
    N_d,id_nd,N_nd,field_names); 
plugin.fields = field_names;            % Create plugin 

options.print = 'yes';                  % Print output to screen (figures)
options.parallel = 'no'; %'yes';               % Multi-core: No IO writing
options.modout = 'yes';                 % Return hmodel simulations
options.save = 'yes';                   % Save memory during trial

% Run the AMALGAM code and obtain non-dominated solution set
[X,F,output,Z] = AMALGAM(AMALGAMPar,Func_name,Par_info,options,plugin);
