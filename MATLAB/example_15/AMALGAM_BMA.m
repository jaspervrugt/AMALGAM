function [FX,G_dot] = AMALGAM_BMA(par,plugin)
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% %
% This function calculates the RMSE, IS, and CRPS of the BMA mixture      %
% distribution                                                            %
%                                                                         %
%  REFERENCES                                                             %
%  Vrugt, J.A. (2024), Distribution-Based Model Evaluation and            %
%      Diagnostics: Elicitability, Propriety, and Scoring Rules for       %
%      Hydrograph Functionals, Water Resources Research, 60,              %
%      e2023WR036710, https://doi.org/10.1029/2023WR036710                %
%  Vrugt, J.A., M.P. Clark, C.G.H. Diks, Q. Duan, and B. A. Robinson      %
%      (2006), Multi-objective calibration of forecast ensembles using    % 
%      Bayesian model averaging, Geophysical Research Letters, 33,        %
%      L19817                                                             %
%                                                                         %
%  © Written by Jasper A. Vrugt, Jan. 2005                                %
%  Los Alamos National Laboratory                                         %
%  University of California Irvine                                        %
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% %

if nargin < 2
    error('AMALGAM_BMA:TooFewInputs',['Requires at least two input ' ...
        'arguments.']);
end

D = plugin.BMA.D; y = plugin.BMA.y; % Ensemble forecasts & verifying data
[n,K] = size(D);                    % # forecasts, # ensmble mmbers
w = par(1:K);                       % Unpack weights

switch plugin.BMA.VAR   % VARIANCE OPTION 
                        % nxK matrix standard deviation forecasts
    case {'1'} % 1: common constant variance
        sigma = par(K+1) * ones(n,K);                                    
    case {'2'} % 2: individual constant variance
        sigma = bsxfun(@times,par(K+1:2*K)',ones(n,K));                               
    case {'3'} % 3: common non-constant variance
        c = par(K+1); sigma = c * D;                               
    case {'4'} % 4: individual non-constant variance
        c = par(K+1:2*K)'; sigma = bsxfun(@times,c,D);
end
sigma = max(sigma,eps);             % sigma >= 2.22e-16

switch plugin.BMA.PDF   % CONDITIONAL DISTRIBUTION 
                        % nxK matrices A and B
    case 'normal'       % Gaussian with µ = D and standard deviation sigma
        A = D; B = sigma; 
    case 'gamma'    % Gamma with shape A and scale B  
        mu = abs(D); var = sigma.^2; A = mu.^2./var; B = var./mu;
end
Y = repmat(y,1,K);                  % Make K copies verifying data
L = pdf(plugin.BMA.PDF,Y,A,B);      % nxK matrix likelihoods forecasts
lik = L*w + realmin;                % nx1 vector likelihoods BMA model
G_dot = D*w; res = y - G_dot;       % nx1 vectors BMA µ frcast & residuals
G = BMA_rnd(plugin.BMA.PDF,w,A,B,2);% Draw 2 nx1-vectors BMA forecasts at y
                                    % Check CRPS using MODELAVG toolbox

FX(1) = sqrt(1/n*(res'*res));       % F_1: RMSE mean forecast BMA model
FX(2) = -mean(log(lik));            % F_2: Ignorance score BMA model
FX(3) = mean(abs(G{1}-y)) -.5 * ... % F_3: CRPS score BMA model
    mean(abs(G{1}-G{2}));

end
