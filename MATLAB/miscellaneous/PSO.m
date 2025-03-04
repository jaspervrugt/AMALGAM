function [G,PS] = PSO(AMALGAMPar,X,PS,n,id)
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% %
% Generate offspring using Paticle Swarm method                           %
%                                                                         %
%  SYNOPSIS                                                               %
%   [G,PS] = PSO(AMALGAMPar,X,PS,n,id)                                    %
%  where                                                                  %
%                                                                         %
%  © Written by Jasper A. Vrugt, Jan. 2005                                %
%  Los Alamos National Laboratory                                         %
%  University of California Irvine                                        %
%                                                                         %
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% %

% Draw random numbers
r_1 = rand(AMALGAMPar.N,AMALGAMPar.d); 
r_2 = rand(AMALGAMPar.N,AMALGAMPar.d);
% Draw varphi
varphi = AMALGAMPar.varphi(AMALGAMPar.N);
% Compute velocity of swarm
PS.v = bsxfun(@times,PS.v,varphi) + AMALGAMPar.c_1 * r_1 .* ...
    ( PS.p - X ) + AMALGAMPar.c_2 * r_2 .* ( PS.n - X );
% Update particle positions
G = X(id,1:AMALGAMPar.d) + PS.v(id,1:AMALGAMPar.d);
% Rotate
xi = 1 + unifrnd(-1,1,n,1); 
% Now  adapt children by random rotation
G = bsxfun(@times,G,xi);

end
