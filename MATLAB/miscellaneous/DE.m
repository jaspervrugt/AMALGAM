function G = DE(AMALGAMPar,X,n,id)
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% %
% Generate offspring using differential evolution                         %
%                                                                         %
%  SYNOPSIS                                                               %
%   G = DE(AMALGAMPar,X,n,id)                                             %
%  where                                                                  %
%                                                                         %
%  © Written by Jasper A. Vrugt, Jan. 2005                                %
%  Los Alamos National Laboratory                                         %
%  University of California Irvine                                        %
%                                                                         %
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% %

% Compute multiplier F
F = bsxfun(@times,AMALGAMPar.beta_1(n),ones(1,AMALGAMPar.d));
% Compute multiplier K
K = bsxfun(@times,AMALGAMPar.beta_2(n),ones(1,AMALGAMPar.d));
% Now create rr values
[~,rr] = sort(rand(n,AMALGAMPar.N),2);
% Generate children using differential evolution
G = X(id,1:AMALGAMPar.d) + F .* (X(rr(:,1),1:AMALGAMPar.d) ...
        - X(id,1:AMALGAMPar.d)) + K .* (X(rr(:,2),1:AMALGAMPar.d) ...
        - X(rr(:,3),1:AMALGAMPar.d));
% Y = X(rr(:,1),1:AMALGAMPar.d) + ...
%     F .* (X(rr(:,1),1:AMALGAMPar.d) - X(rr(:,2),1:AMALGAMPar.d));

end
