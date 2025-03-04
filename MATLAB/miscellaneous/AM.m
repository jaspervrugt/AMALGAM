function G = AM(AMALGAMPar,X,R,n,id)
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% %
% Generate offspring using adaptive Metropolis algorithm                  %
%                                                                         %
%  SYNOPSIS                                                               %
%   G = AM(AMALGAMPar,X,R,n,id)                                           %
%  where                                                                  %
%                                                                         %
%  © Written by Jasper A. Vrugt, Jan. 2005                                %
%  Los Alamos National Laboratory                                         %
%  University of California Irvine                                        %
%                                                                         %
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% %

% Calculate the covariance structure
R = chol(cov(X(R == 1,1:AMALGAMPar.d)) + 1e-10 * eye(AMALGAMPar.d));

% Create children
G = X(id,1:AMALGAMPar.d) + randn(n,AMALGAMPar.d) * ...
    AMALGAMPar.gamma * R;

end
