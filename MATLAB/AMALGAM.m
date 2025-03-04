function varargout = AMALGAM(AMALGAMPar,Func_name,Par_info,varargin)
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% %
%                                                                         %
%     AAA    MMM    MMM    AAA    LLL      GGGGGGGGG    AAA    MMM    MMM %
%    AA AA   MMM    MMM   AA AA   LLL      GGGGGGGGG   AA AA   MMM    MMM %
%   AAA AAA  MMM    MMM  AAA AAA  LLL      GGG   GGG  AAA AAA  MMM    MMM %
%  AAA   AAA MMMM  MMMM AAA   AAA LLL      GGG   GGG AAA   AAA MMMM  MMMM %
%  AAA   AAA MMMMMMMMMM AAA   AAA LLL      GGGGGGGGG AAA   AAA MMMMMMMMMM %
%  AAAAAAAAA MMMMMMMMMM AAAAAAAAA LLL      GGGGGGGGG AAAAAAAAA MMMMMMMMMM %
%  AAAAAAAAA MMM    MMM AAAAAAAAA LLL            GGG AAAAAAAAA MMM    MMM %
%  AAA   AAA MMM    MMM AAA   AAA LLL            GGG AAA   AAA MMM    MMM %
%  AAA   AAA MMM    MMM AAA   AAA LLLLLLLL       GGG AAA   AAA MMM    MMM %
%  AAA   AAA MMM    MMM AAA   AAA LLLLLLLL GGGGGGGGG AAA   AAA MMM    MMM %
%                                                                         %
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% %
%                                                                         %
% AMALGAM: This general purpose MATLAB code is designed to find parameter %
% values that defines the Pareto trade-off surface corresponding to a     %
% vector of different objective functions. In principle, each Pareto      %
% solution is a different weighting of the objectives used. Therefore,    %
% one could use multiple trials with a single objective optimization      %
% algorithms using diferent values of the weights to find different       %
% Pareto solutions. However, various contributions to the optimization    %
% literature have demonstrated that this approach is rather inefficient.  %
% The AMALGAM code developed herein is designed to find an approximation  %
% of the Pareto solution set within a single optimization run. The        %
% AMALGAM method combines two new concepts, simultaneous multimethod      %
% search, and self-adaptive offspring creation, to ensure a fast,         %
% reliable, and computationally efficient solution to multiobjective      %
% optimization problems. This method is called a multi-algorithm,         %
% genetically adaptive multiobjective, or AMALGAM, method, to evoke the   %
% image of a procedure that blends the attributes of the best available   %
% individual optimization algorithms                                      %
%                                                                         %
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% %
%                                                                         %
%  SYNOPSIS                                                               %
%   [X,FX,output,Z,YX] = AMALGAM(AMALGAMPar,Func_name,Par_info);          %
%   [X,FX,output,Z,YX] = AMALGAM(AMALGAMPar,Func_name,Par_info, ...       %
%       options);                                                         %
%   [X,FX,output,Z,YX] = AMALGAM(AMALGAMPar,Func_name,Par_info, ...       %
%       options,plugin);                                                  %
%   [X,FX,output,Z,YX] = AMALGAM(AMALGAMPar,Func_name,Par_info, ...       %
%       options,plugin,Ftrue);                                            %
%  where                                                                  %
%   AMALGAMPar  [input] Structure with AMALGAM settings/parameters        %
%    .d             Number of parameters                                  %
%    .T             Number of generations                                 %
%    .m             Number of objective functions                         %
%    .N             Population size                                       %
%    .rec_methods   Recombination methods                                 %
%      = 'ga'       Genetic algorithm                                     %
%      = 'am'       Adaptive Metropolis                                   %
%      = 'de'       Differential evolution                                %
%      = 'ps'       Particle swarm optimization                           %
%    .beta1         Scaling factor for DE algorithm       DEF: U[0.6,1]   %
%    .beta2         Scaling factor for DE algorithm       DEF: U[0.2,0.6] %
%    .c1            Social factor, PSO algorithm          DEF: 1.5        %
%    .c2            Cognitive factor, PSO algorithm       DEF: 1.5        %
%    .varphi        Inertia factor, PSO algorithm         DEF: U[0.5,1.0] %
%    .p_CR          Crossover probability, GA algorithm   DEF: 0.9        %
%    .p_M           Mutation probability, GA algorithm    DEF: 1/d        %
%    .eta_C         Crossover index, GA algorithm         DEF: 10         %
%    .eta_M         Mutation index, GA algorithm          DEF: 50         %
%    .gamma         Jumprate AMS algorithm                DEF: 2.38^2/d   %
%    .K             Thinning constant                     DEF: 1          %
%    .p0            Minimum probability crossover methods DEF: 0.05       %
%   Func_name   [input] Name (string) function computes objective functs. %
%   Par_info    [input] Parameter structure: Ranges, initial/prior & bnd  %
%    .names         1xd-cell array with parameter names   DEF: []         %
%    .min           1xd-vector of min parameter values    DEF: -inf(1,d)  %
%    .max           1xd-vector of max parameter values    DEF: inf(1,d)   %
%    .boundhandling Treat the parameter bounds or not?                    %
%      = 'reflect'  Reflection method                                     %
%      = 'bound'    Set to bound                                          %
%      = 'fold'     Folding [Vrugt&Braak: doi:10.5194/hess-15-3701-2011]  %
%      = 'none'     No boundary handling                  DEFault         %
%    .initial       Method to draw initial chain states                   %
%      = 'uniform'  Uniform: U(Par_info.min,Par_info.max)                 %
%      = 'latin'    Latin hypercube: LH(Par_info.min,Par_info.max)        %
%      = 'normal'   Normal:  N(Par_info.mu,Par_info.cov)                  %
%      = 'prior'    User specified prior distribution                     %
%      = 'user'     Initial population taken from Par_info.x0             %
%    .mu            1xd-mean vector: µ if .initial = 'normal'             %
%    .cov           dxd-covariance matrix: Σ if .initial = 'normal'       %
%    .prior         Prior distribution (manual) if .initial = 'prior'     %
%                   Ex 1: Par_info.prior = @(x,a,b) mvnpdf(x,a,b);        %
%                             Par_info.a = [-2 -2];                       %
%                             Par_info.b = eye(2);                        %
%                   Ex 2: Par_info.prior = {'normpdf(x,0,1)',...          %
%                                           'unifpdf(x,-2,2)'}            %
%    .steps         d-vector with # intervals for each parameter          %
%                    → Discrete AMALGAM                                   %
%   options     [input] Structure with computational settings/options     %
%    .parallel      Multi-core computation chains?        DEF: 'yes'      %
%    .IO            If parallel, IO writing model?        DEF: 'no'       %
%    .screen        Print screen output during trial?     DEF: 'no'       %
%    .ranking       ABC epsilon value (scalar/vector)                     %
%     = 'matlab'    Nondominated sorting in MATLAB        DEFault         %
%     = 'C++'       Nondominated sorting in C++ (faster)                  %
%    .density       Density of rank 1,2,3, etc. solutions DEF: 'crowding' %
%     = 'crowding'  Crowding distance:Deb et al: NSGA-II  DEFault         %
%     = 'strength'  Strength Pareto:Zitzler&Thiele: SPEA-2                %
%    .modout        Return model simulations?             DEF: 'no'       %
%    .save          Save AMALGAM output during the run?   DEF: 'no'       %
%    .restart       Restart run? (only with "save")       DEF: 'no'       %
%    .print         Output writing screen (tables/figs)   DEF: 'yes'      %
%   plugin      [input] OPT: 2nd argument of Func_name: Class set by user %
%   Ftrue       [input] OPT: NxM matrix true Pareto front, objectve space %
%   X           [outpt] Nxd matrix final population                       %
%   FX          [outpt] Nxm matrix final objective functions              %
%   output      [outpt] Structure summarizes algorithmic performance      %
%    .p_alg         Selection probability crossover methods               %
%    .IGD           Inverse generational distance [if F_par defined]      %
%    .RunTime   [outpt] CPU time in seconds                               %
%   Z           [outpt] NxTxd matrix of past populations                  %
%   YX          [outpt] Model simulations of Pareto solutions             %
%                                                                         %
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% %
%                                                                         %
%  ALGORITHM HAS BEEN DESCRIBED IN                                        %
%   Vrugt, J.A., Multi-criteria optimization using the AMALGAM software   %
%       package: Theory, concepts, and MATLAB implementation, UCI, 2015   %
%   Vrugt, J.A., B.A. Robinson, and J.M. Hyman (2009), Self-adaptive      %
%       multimethod search for global optimization in real-parameter      %
%       spaces, IEEE Transactions on Evolutionary Computation, 13(2),     %
%       pp. 243-259, https://doi.org/10.1109/TEVC.2008.924428             %
%   Vrugt, J.A., and B.A. Robinson (2007), Improved evolutionary          %
%       optimization from genetically adaptive multimethod search,        %
%       Proceedings of the National Academy of Sciences of the United     %
%       States of America, 104, pp. 708-711,                              %
%       https://doi.org/10.1073/pnas.061047110407                         %
%  MORE INFORMATION IN                                                    %
%   Vrugt, J.A., H.V. Gupta, L.A. Bastidas, W. Bouten, and S. Sorooshian  %
%       (2003), Effective and efficient algorithm for multi-objective     %
%       optimization of hydrologic models, Water Resources Research,      %
%       39(8), art. No. 1214, https://doi.org/10.1029/2002WR001746        %
%   Schoups, G.H., J.W. Hopmans, C.A. Young, J.A. Vrugt, and              %
%       W.W. Wallender, Multi-objective optimization of a regional        %
%       spatially-distributed subsurface water flow model, Journal of     %
%       Hydrology, pp. 20-48, 311(1-4),                                   %
%       https://doi.org/10.1016/j.jhydrol.2005.01.001, 2005.              %
%   Vrugt, J.A., P.H. Stauffer, T. Wöhling, B.A. Robinson, and            %
%       V.V. Vesselinov (2008), Inverse modeling of subsurface flow and   %
%       transport properties: A review with new developments, Vadose      %
%       Zone Journal, 7(2), 843 - 864,                                    %
%       https://doi.org/10.2136/vzj2007.0078                              %
%   Wöhling, T., J.A. Vrugt, and G.F. Barkle (2008), Comparison of three  %
%       multiobjective optimization algorithms for inverse modeling of    %
%       vadose zone hydraulic properties, Soil Science Society of America %
%       Journal, 72, 305 - 319, https://doi.org/10.2136/sssaj2007.0176    %
%   Wöhling, T., and J.A. Vrugt (2008), Combining multi-objective         %
%       optimization and Bayesian model averaging to calibrate forecast   %
%       ensembles of soil hydraulic models, Water Resources Research, 44, %
%       W12432, https://doi.org/10.1029/2008WR007154                      %
%                                                                         %
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  %
%                                                                         %
%  BUILT-IN CASE STUDIES                                                  %
%   Example 1   Multivariate normal benchmark study                       %
%   Example 1   ZDT1: test function                                       %
%   Example 2   ZDT2: test function                                       %
%   Example 3   ZDT3: test function                                       %
%   Example 4   ZDT4: test function                                       %
%   Example 5   ZDT6: test function                                       %
%   Example 6   ZDT6: test function, discrete parameter space             %
%   Example 7   DTLZ1: test function, 3 objectives                        %
%   Example 8   DTLZ2: test function, 3 objectives                        %
%   Example 9   DTLZ3: test function, 3 objectives                        %
%   Example 11  Real-world example using rainfall-discharge modeling      %
%   Example 12  Watershed modeling using driven & nondriven hydrograph    %
%   Example 13  Bayesian model averaging: RMSE, IS and CRPS               %
%   Example 14  Multi-criteria BMA training temperature ensemble          %
%   Example 15  Multi-criteria BMA training sea-level pressure ensemble   %
%                                                                         %
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% %
%                                                                         %
%  © Written by Jasper A. Vrugt, Jan. 2005                                %
%  Los Alamos National Laboratory                                         %
%  University of California Irvine                                        %
%  Version 0.5    June 2006                                               %
%  Version 1.0    January 2009    Cleaning & implemented test problems    %
%  Version 1.1    January 2010    Variable population size                %
%  Version 1.2    August 2010     Sampling from prior distribution        %
%  Version 1.3    May 2014        Cleaning and improved speed ranking     %
%  Version 1.4    Januari 2014    Parallellization (done if CPU > 1)      %
%  Version 2.0    June 2024       Recoding - and other updates            %
%                                                                         %
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% %

% Check the input arguments
if nargin < 3, error(['AMALGAM ERROR: Main function of AMALGAM ' ...
        'requires at least three input arguments.']); end
options = struct('restart','no'); plugin = struct; Ftrue = []; 
if nargin == 4, options = varargin{1}; end
if nargin == 5, [options,plugin] = varargin{1:2}; end
if nargin == 6, [options,plugin,Ftrue] = varargin{1:3}; end
if ~isfield(options,'restart'), options.restart = 'no'; end
% <><><><><><><><><><><>< END CHECK INPUT ARGS ><><><><><><><><><><><><><>

% <><><><><><><><><><><><><> INITIALIZATION <><><><><><><><><><><><><><><>
switch options.restart
    case 'no' % New trial: No restart
        [AMALGAMPar,Par_info,options] = ...         % Check, stop if wrong
            AMALGAM_check(Func_name,...         
            AMALGAMPar,Par_info,options,Ftrue); 
        [AMALGAMPar,Par_info,options, ...           % Check settings user
            T_start] = AMALGAM_setup(...   
            AMALGAMPar,Par_info,options);
        [AMALGAMPar,f_handle] = ...                 % Computational setup 
            AMALGAM_calc_setup(AMALGAMPar, ...       
            Func_name,options,plugin);
        [AMALGAMPar,Par_info,X,p_rm,PS,Z, ...       % Initialize variables
            output] = AMALGAM_initialize( ...              
            AMALGAMPar,Par_info);
        [FX,YX] = AMALGAM_calc_FX(X,AMALGAMPar, ... % Compute OFs parents
            f_handle,options);
        [RX,dX,FX_min] = AMALGAM_rank(FX, ...       % Rank parents & distance
            options);     
        Z(1:AMALGAMPar.N,1:AMALGAMPar.d + ...       % Add X to Z
            AMALGAMPar.m) = [X FX];      
        if any(strcmp(AMALGAMPar.rec_methods,'ps')) % Update PS structure
            PS = Update_PS(AMALGAMPar, ...             
                    Z(1:AMALGAMPar.N,1:AMALGAMPar.d ...
                    + AMALGAMPar.m),PS,FX_min);
        end     
        if ~isempty(Ftrue)                          % Compute hypervolume 
            output.IGD(1,1:2) = [1 IGD(Ftrue',FX')]; 
        end 
    case 'yes' % Restart run
        [AMALGAMPar,Par_info,options,PS, ...        % Load restart variabls
            f_handle,X,output,Z,FX,YX,...
            FX_min,RX,dX,T_start] = ...
            AMALGAM_restart('AMALGAM.mat');  
end
cnt = 0; ct = 1; t0 = cputime; [fg,ax0] = deal([]); % Start timer

% <><><><><><><><><><><><> END INITIALIZATION <><><><><><><><><><><><><><>

% <><><><><><><><><><><><><><> DYNAMIC PART <><><><><><><><><><><><><><><>
for t = T_start : AMALGAMPar.T                          % Now iterate
    [id,id_rm] = AMALGAM_distribution(AMALGAMPar,p_rm); % Create distribution rec. methods
    G = AMALGAM_children(AMALGAMPar,Par_info,X,FX, ...  % Generate offspring
        RX,dX,PS,id_rm);
    [FG,YG] = AMALGAM_calc_FX(G,AMALGAMPar,f_handle,... % Objective functions children 
        options);
    FXG_min = min(FX_min,min(FG));                      % Minimum objective functions
    [X,FX,RX,dX,id_N,id] = ...                          % New population
        AMALGAM_population(AMALGAMPar,options,X, ...    
        G,FX,FG,id);
    p_rm = AMALGAM_load(AMALGAMPar,p_rm,id);            % New probability rec. methods
    if strcmp(options.modout,'yes')                     % Update simulations current generation
        simtot = [YX;YG]; YX = simtot(id_N,:);
    end
    if mod(t,AMALGAMPar.K) == 0                         
        Z(ct*AMALGAMPar.N+1:(ct+1)*AMALGAMPar.N, ...    % Append to archive 
            1:AMALGAMPar.d+AMALGAMPar.m) = [X FX]; 
            ct = ct + 1;
        if any(strcmp(AMALGAMPar.rec_methods,'ps'))     % Update PS structure
            PS = Update_PS(AMALGAMPar, ...             
                    Z(1:ct*AMALGAMPar.N,1:AMALGAMPar.d ...
                    + AMALGAMPar.m),PS,FXG_min);
        end
    end
    output.p_rm(t,1:AMALGAMPar.q+1) = [t p_rm];         % Store prob. rec. methods
    if ~isempty(Ftrue)                                  % Compute hypervolume
        output.IGD(t,1:2) = [t IGD(Ftrue',FX')];
    end
    if strcmp(options.save,'yes')                       % Save ouput?
        save AMALGAM.mat 
    end            
    if strcmp(options.screen,'yes')                     % Plot to screen? 
        [fg,ax0] = plot_screen(fg,ax0, ...
            AMALGAMPar,FX,t,Ftrue,output); 
    end
    if (t > 1)                                          % Print progress
        fprintf(1, repmat('\b',1,cnt));
        cnt = fprintf(['AMALGAM calculating, ' ...
            '%% done: %3.2f'],100 *(t/AMALGAMPar.T));
    end
end
% <><><><><><><><><><><><> END DYNAMIC PART <><><><><><><><><><><><><><><>
fprintf(1,'\n');                                        % Skip a line

% <><><><><><><><><><><><>< POSTPROCESSING ><><><><><><><><><><><><><><><>
output.RunTime = cputime - t0;                          % Determine total run time
AMALGAM_end(AMALGAMPar,options);                        % Close workers
Z = Z(1:ct*AMALGAMPar.N,1:AMALGAMPar.d+AMALGAMPar.m);   % Finalize achive Z
if strcmp(options.print,'yes')                          % Plot figures
    cnt = 0;
    for t = 2:3                                         % Print progress
        if t == 2, print_name = '........'; end
        if t == 3, print_name = '........ done'; end
        if ( t > 1 )
            fprintf(1, repmat('\b',1,cnt));             % Delete line before
            cnt = fprintf(['AMALGAM ' ...
                'postprocessor, please wait %s'], ...
                print_name);
            if t == 2
                AMALGAM_postproc(AMALGAMPar, ...        % Postprocess results
                    Par_info,options,X,FX, ...
                    output,Ftrue,YX,Z);
            end
        end
    end
end
fprintf(1,'\n');                                        % Skip a line
varargout(1:5) = { X FX output Z YX };                  % Return arguments
% <><><><><><><><><><><>< END POSTPROCESSING <><><><><><><><><><><><><><><

end
