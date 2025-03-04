function FX_true = Compute_FX_true(AMALGAMPar,Func_name,problem,M,plugin)
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% %
%Compute_FX_true: Computes actual Pareto optimal solutions                %
% Note: For ZDT1, ZDT2, ZDT3, ZDT4 and ZDT6                               %
% Must revisit code for PNAS paper, 2007 to obtain Pareto solution set    %
% for other problems                                                      %
%                                                                         %
%  SYNOPSIS                                                               %
%   FX_true = Compute_FX_true(AMALGAMPar,Func_name,Par_info, ...          %
%       options,plugin,M,problem)                                         %
%  where                                                                  %
%   AMALGAMPar  [input] Structure with AMALGAM settings/parameters        %
%   Func_name   [input] Name (string) function computes objective functs. %
%   problem     [input] Name (string) of benchmark problem                %
%       'ZDT1'      ZDT1 of Zitzler, Deb & Thiele 2000                    %
%       'ZDT2'      ZDT2 of Zitzler, Deb & Thiele 2000                    %
%       'ZDT3'      ZDT3 of Zitzler, Deb & Thiele 2000                    %
%       'ZDT4'      ZDT4 of Zitzler, Deb & Thiele 2000                    %
%       'ZDT6'      ZDT6 of Zitzler, Deb & Thiele 2000                    %
%       'DTLZ1'     DTLZ1 of Deb et al. 2001                              %
%   M           [input] OPT: # points true Pareto front  DEF: 500         %
%   plugin      [input] OPT: 2nd argument of Func_name: Class set by user %
%   FX_true     [outpt] True Pareto front (approximation) objective space % 
%                                                                         %
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% %
%                                                                         %
%  © Written by Jasper A. Vrugt, Jan. 2005                                %
%  Los Alamos National Laboratory                                         %
%  University of California Irvine                                        %
%                                                                         %
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% %

options = struct('restart','no','parallel','no','modout','no');
if nargin < 5, plugin = struct; end
if nargin < 4, M = 500; end
if nargin < 3, error(['Compute_FX_true:Requires at least three ' ...
        'input arguments']); end
problem = upper(problem);                       % Capital letters problem

[AMALGAMPar,f_handle] = ...                     % Computational setup 
    AMALGAM_calc_setup(AMALGAMPar, ...       
    Func_name,options,plugin);
d = AMALGAMPar.d; N_max = 50000;                % # points to sample

switch char(upper(problem))
    case {'ZDT1','ZDT2','ZDT3','ZDT4','ZDT6'}   % Zitzler's functions
        FX1_min = zeros(1,AMALGAMPar.d); 
        FX2_min = [1 zeros(1,AMALGAMPar.d-1)];
    case {'DTLZ1','DTLZ2','DTLZ3'}              % Deb's functions
        FX1_min = [0 0 .5 * ones(1,AMALGAMPar.d-2)]; 
        FX2_min = [1 1 .5 * ones(1,AMALGAMPar.d-2)];
    case 'hmodel'
        FX1_min = []; FX2_min = [];
    otherwise
        error('Compute_FX_true:Do not know this function')
end

if ~isempty(FX1_min)
    T = size(FX1_min,1); N_tot = [0 nan(1,T)];  
    for j = 1:T
        N_tot(j+1) = round(N_max/T);            % # samples?
    end
    X_all = nan(sum(N_tot),AMALGAMPar.d);       % Initialize samples
    for j = 1:T
        X = LH_sampling(FX1_min(j,:), ...       % Latin hypercube sampling
                FX2_min(j,:),N_tot(j+1));
        X(1,1:d) = FX1_min(j,1:d);              % Single objective optimum
        X(2,1:d) = FX2_min(j,1:d);              % Single objective optimum
        X_all(N_tot(j)+1:N_tot(j+1),1:d) = X;   % Store samples
    end
    FX_all = AMALGAM_calc_FX(X_all, ...         % Compute objective values
        AMALGAMPar,f_handle,options);
    FX_all = sortrows(FX_all,1);                % Sort FX_all ascndng order
    FX_true = getF_true(FX_all,AMALGAMPar,1,M+1); % Uniform Pareto approx.
    % Always ok except when dealing with Fonseca & Fleming function
end

end
%<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>

%<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>
% Secondary functions
%<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>

% 1: getF_true
function FX_true = getF_true(FX,AMALGAMPar,N_f,M)
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% %
%getF_true: Generates 202 uniformly spaced Pareto solutions
%
% Written by Jasper A. Vrugts
% Los Alamos National Laboratory
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% %

m = AMALGAMPar.m;                       % Extract # objectives
FX = sortrows(FX,1);                    % First sort FX - first objective
id = find(diff(FX(:,1)) > 0);           % Find dissimilar FX points
FX = [FX(1,1:m) ; ...                   % Keep only dissimilar FX points
    FX(id+1,1:m)];
N_FX = size(FX,1);                      % # dissimilar points
Ed = nan(N_FX,1); Ed(1) = 0;            % Eucld distance adjacent FX points
for i = 2:N_FX
    Ed(i) = norm(FX(i,1:m) - FX(i-1,1:m),2);
end
[Ed_s,ii_s] = sort(-Ed);                % Sort Ed in decreasing distance
if N_f == 1
    id = [0 N_FX]'; Ed_cs = cumsum(Ed);
else
    id = [0 , sort(ii_s(1:N_f-1))' , ...
        N_FX + 1]; 
    Ed_cs = cumsum(-Ed_s(N_f:N_FX));
end
mEd = Ed_cs(end)/M; %(M-1);             % Mean distance between points
FX_true = nan(M-2,m);                   % Initialize Pareto approximation
nq = [0 nan(1,N_f)];                    % Recursive addition to F_true
for i = 1:N_f
    Ed_f = cumsum(Ed(id(i)+1:id(i+1)-1,1));     % Define Ed of front
    Ed_cs = Ed_f(end);
    FX_f = FX(id(i)+1:id(i+1)-1,1:m);           % Define FX of front
    nf = size(FX_f,1);                          % # points?
    if mEd > Ed_cs
        F_pareto = FX_f(1,:);                   % For interpolation
    else
        Ed_int = mEd : mEd : Ed_cs; 
        for z = 1:m
            F_pareto(:,z) = interp1(Ed_f,FX_f(:,z),Ed_int);   % Linear interpolation
        end
    end
    n_p = size(F_pareto,1);                     % # points F_pareto?
    F_pareto(1,1:m) = FX_f(1,1:m);              % First point F_pareto
    F_pareto(n_p,1:m) = FX_f(nf,1:m);           % Last point F_pareto
    q = find(~isnan(F_pareto(:,1)));            % Outside bound values
    nq(i+1) = numel(q); nq = cumsum(nq);        % # points & cumulative 
    FX_true(nq(i)+1:nq(i+1),1:m) = ...          % Now add to F_true
        F_pareto(q,1:m);      
end

end
%<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>
