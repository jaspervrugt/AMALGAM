function [Xn,FXn,RXn,dXn,id_N,id] = AMALGAM_population(AMALGAMPar, ...
    options,X,G,FX,FG,id)
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
% Function selects new population based on current offspring and parents  % 
%                                                                         %
%  SYNOPSIS                                                               %
%   [Xn,FXn,RXn,dXn,id_N,id] = AMALGAM_population(AMALGAMPar, ...         %
%       options,X,G,FX,FG,id)                                             %
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

Q = [X ; G]; FQ = [FX ; FG];        % Append children to parents
[RQ,dQ] = AMALGAM_rank(FQ,options); % Rank jointly & crowding distance
I_alg = [zeros(AMALGAMPar.N,1);id]; % Indices of recombination methods

% Selection: Based on rank
RQ_max = max(RQ);                   % Maximum rank of Q
n_rnk = nan(RQ_max,1);              % # points with given rank
for r = 1:RQ_max
    id_r = find(RQ == r);           % id stores index each rank
    n_rnk(r) = numel(id_r);         % # points rank "r"
end
tot_rnk = cumsum(n_rnk);            % Which ranks automatic new generation
idx = find(tot_rnk<=AMALGAMPar.N);  % Which ranks <= spots to fill

if ~isempty(idx)
    r_max = idx(end);               % Maximum rank
    id_R = find(RQ <= r_max);       % Select all points rank = 1 to r_max
    n_lft = AMALGAMPar.N ...        % # remaining spots to fill
        -  tot_rnk(r_max);          
    id_lft = find(RQ == r_max + 1); % Go to next rank
    [~,ii] = sortrows(-dQ(id_lft)); % Sort crowding distance descending order
    id_C = id_lft(ii(1:n_lft));     % Pick largest n_lft distances
    id_N = [id_R ; id_C];           % Combine id_Rank and id_Crowd
else
    id_R = find(RQ == 1);           % Pick rank 1 only but crowding distance
    [~,ii] = sortrows(-dQ(id_R));   % Sort crowding distance descending order
    id_N = id_R(ii(1:AMALGAMPar.N));% Combine id_Rank and id_Crowd
end

% Now extract new population, Xn(ew)
Xn = Q(id_N,1:AMALGAMPar.d);        % Xn
FXn = FQ(id_N,1:AMALGAMPar.m);      % FXn
RXn = RQ(id_N); dXn = dQ(id_N);     % Rank & crowding distnce
id = I_alg(id_N);                   % Nx1 vector indcs recmbntion methds
                                    % 0, ... , AMALGAMPar.q

end
