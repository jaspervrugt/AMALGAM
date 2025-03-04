function G = BMA_rnd(PDF,w,A,B,P)
% This function samples trajectories of BMA forecast distribution

if nargin < 5
    error('BMA_rnd:TooFewInputs',['Requires at least five ' ...
        'input arguments.']);
end

[n,K] = size(A);                        % # forecasts and ensemble members  
G = cell(1,K);                          % Cell array for components
for p = 1:P                             % N trajectories BMA forcst dstrbtion
    id = randsample(1:K,n,'true',w);    % Draw nx from 1:K using weights
    for k = 1:K
        T = find(id == k);              % kth component for indices of T
        G{p}(T,1) = ...                 % Sample from kth mixture component
            random(PDF,A(T,k),B(T,k));  
    end
end