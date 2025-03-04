function rank = FNS(Q)
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% %
% Fast nondominated sorting for Nxm matrix of objective function values   %
%  SYNOPSIS: rank = FNS(Q)                                                %
%   where                                                                 %
%    Q         [input] Mxm matrix objective functions parents & children  %
%    rank      [outpt] Mx1 vector of Pareto ranks N individual points     %
%                                                                         %
%  © Written by Jasper A. Vrugt, Jan. 2005                                %
%  Los Alamos National Laboratory                                         %
%  University of California Irvine                                        %
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% %
[M,m] = size(Q);    % Number of individuals and objective functions
Fr = {[]};          % Pareto-optimal fronts
rank = zeros(M,1);  % Initialize vector with ranks
Sp = cell(M,1);     % Individuals a particular individual dominates
np = zeros(M,1);    % # individuals by which a particular point is dominated
for j = 1 : M       % Loop over all elements of R
    idx = find( (sum( bsxfun(@le,Q(j,1:m),Q) , 2 ) == m ) & ...
        ( sum( bsxfun(@lt,Q(j,1:m),Q) , 2 ) > 0 ) );     % Which of Q, j dominates
    if numel(idx), Sp{j} = idx; end                      % Store points in Sp
    idx = find( (sum(bsxfun(@le,Q,Q(j,1:m)) , 2) == m ) & ...
        ( sum(bsxfun(@lt,Q,Q(j,1:m)) , 2 ) > 0 ) );      % Which of Q dominate j
    np(j) = np(j) + numel(idx);                          % Store # points in np
    if np(j) == 0                                        % j is in first front
        Fr{1}(end+1) = j; 
    end
end
i = 1;
while ~isempty(Fr{i})                                    % Continue if not all ranked
    NextFr = [];                                         % Next front is empty
    for j = 1 : numel(Fr{i})                             % For each member j of Fr{i}
        q = Sp{Fr{i}(j)};                                % Modify each member set Sp
        np(q) = np(q) - 1;                               % Decrement np(q) by one
        id = find(np(q) == 0); NextFr = [NextFr;q(id)];  % q member next front
    end
    i = i + 1;                                           % Go to next front
    Fr{i} = NextFr;                                      % Current members of NextFr
end
for j = 1:i, rank(Fr{j}) = j; end                        % Extract M-vector of ranks 

end