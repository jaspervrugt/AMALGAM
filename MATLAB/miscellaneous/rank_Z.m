function id = rank_Z(Z,AMALGAMPar,options)
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% %
%rank_Z: Function which ranks all past populations iteratively
% Written by Jasper A. Vrugt
% University of California Irvine
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% %

[N,d2] = size(Z);               % Determine size of Z
d = d2 - AMALGAMPar.m;          % # parameters
m = AMALGAMPar.m;               % # objective functions
Nrank = 1e3;                    % Maximum # points to rank at once

switch N <= Nrank
    case 1      % Rank Z at once
        R = AMALGAM_rank(Z,options); id = find(R == 1);
    otherwise   % Take segments of Z
        q = floor(N/Nrank); R1 = zeros(N,1) > 0; 
        for ii = 1:q
            id_s = (ii - 1) * Nrank + 1; id_e = ii * Nrank;
            R = AMALGAM_rank(Z(id_s:id_e,d+1:d+m),options);
            R1(id_s:id_e,1) = R == 1;
        end
        if id_e < N
            R = AMALGAM_rank(Z(id_e+1:N,d+1:d+m),options);
            R1(id_e+1:N,1) = R == 1;
        end
        % Locate indices of rank 1 solutions
        id_1 = find(R1 == 1);
        % Now rank all rank 1 solutions of different segments of Z at once
        R = AMALGAM_rank(Z(id_1,d+1:d+m),options);
        % Which are now the rank 1 solutions of Z?
        id_1F = R==1; id = id_1(id_1F);
end