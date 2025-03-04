function [G,PS] = AMALGAM_children(AMALGAMPar,Par_info,X,FX, ...
    RX,dX,PS,id_rm,plugin)
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
% This function uses different recombination methods to generate children %
%                                                                         %
%  SYNOPSIS                                                               %
%   G = AMALGAM_children(AMALGAMPar,Par_info,X,FX,RX,dX,PS,id_rm,plugin)  %
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

% First initialize offspring
G = nan(AMALGAMPar.N,AMALGAMPar.d);

% <><><><><><><><><><><> NSGA-II GENETIC ALGORITHM <><><><><><><><><><><><>
run_GA = find(strcmpi(AMALGAMPar.rec_methods,'ga'));
if run_GA
    id_GA = id_rm{run_GA};                      % Indices of population?
    G(id_GA,1:AMALGAMPar.d) = ...               % Children genetic algorithm
        NSGA(AMALGAMPar,Par_info,X,FX,id_GA,dX);
end

% <><><><><><><><><><> ADAPTIVE METROPOLIS ALGORITHM <><><><><><><><><><><>
run_AM = find(strcmpi(AMALGAMPar.rec_methods,'am'));
if run_AM
    n_AM = numel(id_rm{run_AM});                % # children AM must create
    id_AM = id_rm{run_AM};                      % Indices of population
    G(id_AM,1:AMALGAMPar.d) = ...               % Children adptve Metrpolis
        AM(AMALGAMPar,X,RX,n_AM,id_AM);
end

% <><><><><><><><><><><><> PARTICLE SWARM METHOD <><><><><><><><><><><><><>
run_PS = find(strcmpi(AMALGAMPar.rec_methods,'ps'));
if run_PS
    n_PS = numel(id_rm{run_PS});                % # children PS must create
    id_PS = id_rm{run_PS};                      % Indices of population   
    [G(id_PS,1:AMALGAMPar.d),PS] = ...          % Children Particle Swarm
        PSO(AMALGAMPar,X,PS,n_PS,id_PS);
%    PS.v = V;                                   % Update velocity
end

% <><><><><><><><><><><><> DIFFERENTIAL EVOLUTION <><><><><><><><><><><><><
run_DE = find(strcmpi(AMALGAMPar.rec_methods,'de'));
if run_DE
    n_DE = numel(id_rm{run_DE});                % # children DE must create
    id_DE = id_rm{run_DE};                      % Indices of population    
    G(id_DE,1:AMALGAMPar.d) = ...               % Children Diff. Evolution
        DE(AMALGAMPar,X,n_DE,id_DE);
end

% Boundary handling
if isfield(Par_info,'boundhandling')
    [G,v] = Boundary_handling(G,Par_info);      % v is a Nx1 vector 
else                                            % 0 in bound, 1 otherwise
    v = ones(AMALGAMPar.N,1) < 0;
end
% BMA model training
if isfield(Par_info,'unit_simplex')   
    w_sum = sum(G(1:AMALGAMPar.N,1:plugin.BMA.K),2);    % Sum of BMA weights
    G(:,1:plugin.BMA.K) = bsxfun(@rdivide, ...         	% Weights unit Simplex
        G(:,1:plugin.BMA.K),w_sum);    
end
% Transform to discrete space
if isfield(Par_info,'steps')                    
    G = Discrete_space(G,Par_info);             % Discrete space
end

end
