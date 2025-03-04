function [AMALGAMPar,Par_info,X,p_rm,PS,Z,output] = ...
    AMALGAM_initialize(AMALGAMPar,Par_info,plugin)
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
%  SYNOPSIS                                                               %
%   [AMALGAMPar,Par_info,X,p_rm,PS,Z,output] = AMALGAM_initialize( ...    %
%       AMALGAMPar,Par_info )                                             %
%  where                                                                  %
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
%                                                                         %
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% %
%                                                                         %
%  © Written by Jasper A. Vrugt, Jan. 2005                                %
%  Los Alamos National Laboratory                                         %
%  University of California Irvine                                        %
%                                                                         %
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% %

AMALGAMPar.q = numel(AMALGAMPar.rec_methods);   % # recombination methods
p_rm = 1/AMALGAMPar.q * ones(1,AMALGAMPar.q);   % Selct prob. rec. methods
AMALGAMPar.p0 = max(2/AMALGAMPar.N,5e-2);       % Min prob. rec. methods

output.p_rm = nan(AMALGAMPar.T,AMALGAMPar.q+1); % Init. matrix with p_alg 
output.p_rm(1,1:AMALGAMPar.q+1) = [1 p_rm];     % Store p_alg rec. methods
output.IGD = nan(AMALGAMPar.T,2);               % Init. matrix inv. gen. d 
                                                % if Pareto front known
if any(strcmp(AMALGAMPar.rec_methods,'ps'))     % Particle Swarm used
    % vPSO = rand(AMALGAMPar.N,AMALGAMPar.d);
    PS.v = 1/5 * bsxfun(@times,unifrnd(-1,1, ...
        AMALGAMPar.N,AMALGAMPar.d),Par_info.max - Par_info.min );
else
    PS = struct; 
end

% Initialize initial population
X = nan(AMALGAMPar.N,AMALGAMPar.d);
% Fill spots of initial population
switch Par_info.initial
    case 'uniform'  % Uniform random sampling
        X = Par_info.min + rand(AMALGAMPar.N,AMALGAMPar.d) .* ...
                (Par_info.max - Par_info.min);
    case 'latin'    % Latin hypercube sampling
        X = LH_sampling(Par_info.min,Par_info.max,AMALGAMPar.N);
    case 'normal'   % Multi-normal distribution
        X = repmat(Par_info.mu,AMALGAMPar.N,1) + randn(AMALGAMPar.N, ...
            AMALGAMPar.d) * chol(Par_info.cov);
    case 'prior'    % User-defined prior distribution
        switch Par_info.u
            case 'yes'
                % Univariate prior: Draw one parameter at a time
                for qq = 1:AMALGAMPar.d
                    for zz = 1:AMALGAMPar.N
                        X(zz,qq) = Par_info.prior_rnd{qq}(1);
                    end
                end
            case 'no'
                % Multivariate prior: Draw all parameters at once
                for zz = 1:AMALGAMPar.N
                    X(zz,1:AMALGAMPar.d) = Par_info.prior_rnd(1);
                end
        end
    case 'user' % Initial population specified by user
        for zz = 1:AMALGAMPar.N
            X(zz,1:AMALGAMPar.d) = Par_info.x0(zz,1:AMALGAMPar.d);
        end    
    otherwise   % Display error warning
        error('AMALGAM_initialize:Unknown initial sampling');
end

% Boundary handling
if isfield(Par_info,'boundhandling')
    [X,v] = Boundary_handling(X,Par_info);          %#ok v is a Nx1 vector 
else                                                % 0 in bound, 1 otherwise
    v = ones(AMALGAMPar.N,1) < 0;                   %#ok, REVISE&USE LATER
end

% BMA model training
if isfield(Par_info,'unit_simplex')   
    wght_sum = sum(X(1:AMALGAMPar.N,1:plugin.BMA.K),2);    % Sum of BMA weights
    X(:,1:plugin.BMA.K) = bsxfun(@rdivide, ...             % Weights unit Simplex
        X(:,1:plugin.BMA.K),wght_sum);    
end

% Transform to discrete space
if isfield(Par_info,'steps')                        % Transform discrete space
    X = Discrete_space(X,Par_info);
end
Z = nan(AMALGAMPar.T * AMALGAMPar.N / ...           % Allocate memory Z    
        AMALGAMPar.K, AMALGAMPar.d + AMALGAMPar.m);

% Print to screen
fprintf('\n'); fprintf('\n'); 
fprintf('  -----------------------------------------------------------------------------  \n');
fprintf('        AAA    MMM    MMM    AAA    LLL      GGGGGGGGG    AAA    MMM    MMM      \n');
fprintf('       AA AA   MMM    MMM   AA AA   LLL      GGGGGGGGG   AA AA   MMM    MMM      \n');
fprintf('      AAA AAA  MMM    MMM  AAA AAA  LLL      GGG   GGG  AAA AAA  MMM    MMM      \n');
fprintf('     AAA   AAA MMMM  MMMM AAA   AAA LLL      GGG   GGG AAA   AAA MMMM  MMMM      \n');
fprintf('     AAA   AAA MMMMMMMMMM AAA   AAA LLL      GGGGGGGGG AAA   AAA MMMMMMMMMM      \n');
fprintf('     AAAAAAAAA MMMMMMMMMM AAAAAAAAA LLL      GGGGGGGGG AAAAAAAAA MMMMMMMMMM      \n');
fprintf('     AAAAAAAAA MMM    MMM AAAAAAAAA LLL            GGG AAAAAAAAA MMM    MMM      \n');
fprintf('     AAA   AAA MMM    MMM AAA   AAA LLL            GGG AAA   AAA MMM    MMM      \n');
fprintf('     AAA   AAA MMM    MMM AAA   AAA LLLLLLLL       GGG AAA   AAA MMM    MMM      \n');
fprintf('     AAA   AAA MMM    MMM AAA   AAA LLLLLLLL GGGGGGGGG AAA   AAA MMM    MMM      \n');
fprintf('  -----------------------------------------------------------------------------  \n');
fprintf('  © Jasper A. Vrugt, University of California Irvine \n');
fprintf('\n'); 

end