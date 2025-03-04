function [RQ,dQ,FQ_min] = AMALGAM_rank(FQ,options)
% % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% %
% %                                                                         %
% %     AAA    MMM    MMM    AAA    LLL      GGGGGGGGG    AAA    MMM    MMM %
% %    AA AA   MMM    MMM   AA AA   LLL      GGGGGGGGG   AA AA   MMM    MMM %
% %   AAA AAA  MMM    MMM  AAA AAA  LLL      GGG   GGG  AAA AAA  MMM    MMM %
% %  AAA   AAA MMMM  MMMM AAA   AAA LLL      GGG   GGG AAA   AAA MMMM  MMMM %
% %  AAA   AAA MMMMMMMMMM AAA   AAA LLL      GGGGGGGGG AAA   AAA MMMMMMMMMM %
% %  AAAAAAAAA MMMMMMMMMM AAAAAAAAA LLL      GGGGGGGGG AAAAAAAAA MMMMMMMMMM %
% %  AAAAAAAAA MMM    MMM AAAAAAAAA LLL            GGG AAAAAAAAA MMM    MMM %
% %  AAA   AAA MMM    MMM AAA   AAA LLL            GGG AAA   AAA MMM    MMM %
% %  AAA   AAA MMM    MMM AAA   AAA LLLLLLLL       GGG AAA   AAA MMM    MMM %
% %  AAA   AAA MMM    MMM AAA   AAA LLLLLLLL GGGGGGGGG AAA   AAA MMM    MMM %
% %                                                                         %
% % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% %
% %                                                                         %
% % Pareto ranking and crowding distance computation of solutions matrix FQ %
% % The "quality" of each individual of FQ is summarized by two metrics     %
% %   1. Their rank ( = number of their Pareto-optimal front )              %
% %   2. Their density ( = closeness of neighboring solutions )             %
% % The "best" solution has rank = 1 and infinite distance                  %
% %                                                                         %
% %  SYNOPSIS                                                               %
% %   [R,d] = AMALGAM_rank(FQ)                                              %
% %   [R,d] = AMALGAM_rank(FQ,options)                                      %
% %  where                                                                  %
% %   FQ          [input] 2Nxm matrix of objective function values          %
% %   options     [input] OPT: AMALGAM input structure algorithmic settings %
% %    .density       Density of rank 1,2,3, etc. solutions DEF: 'crowding' %
% %     = 'crowding'  Crowding distance:Deb et al: NSGA-II  DEFault         %
% %     = 'strength'  Strength Pareto:Zitzler&Thiele: SPEA-2                %
% %    .ranking       ABC epsilon value (scalar/vector)                     %
% %     = 'matlab'    Nondominated sorting in MATLAB        DEFault         %
% %     = 'C++'       Nondominated sorting in C++ (faster)                  %
% %   R           [outpt] 2Nx1 vector of rank of points FQ                  %
% %   d           [outpt] 2Nx1-vector of distance of points FQ              %
% %                                                                         %
% % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% %
% %                                                                         %
% %  ALGORITHM HAS BEEN DESCRIBED IN                                        %
% %   Vrugt, J.A., Multi-criteria optimization using the AMALGAM software   %
% %       package: Theory, concepts, and MATLAB implementation, UCI, 2015   %
% %   Vrugt, J.A., B.A. Robinson, and J.M. Hyman (2009), Self-adaptive      %
% %       multimethod search for global optimization in real-parameter      %
% %       spaces, IEEE Transactions on Evolutionary Computation, 13(2),     %
% %       pp. 243-259, https://doi.org/10.1109/TEVC.2008.924428             %
% %   Vrugt, J.A., and B.A. Robinson (2007), Improved evolutionary          %
% %       optimization from genetically adaptive multimethod search,        %
% %       Proceedings of the National Academy of Sciences of the United     %
% %       States of America, 104, pp. 708-711,                              %
% %       https://doi.org/10.1073/pnas.061047110407                         %
% %                                                                         %
% % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% %
% %                                                                         %
% %  © Written by Jasper A. Vrugt, Jan. 2005                                %
% %  Los Alamos National Laboratory                                         %
% %  University of California Irvine                                        %
% %                                                                         %
% % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% %
% 
if nargin < 2
    % We resort to default settings
    warning(['AMALGAM_rank WARNING: Density operator and ranking ' ...
        'language not specified -> crowding distance of ' ...
        'Deb et al. (1990) in matlab used as default']);
    % Define density calculation
    options.density = 'crowding'; options.ranking = 'matlab';
end

[M,m] = size(FQ);        % # solutions Q (M = 2N) and objective functions
RQ = nan(M,1);           % Initialize vector with ranks
FQ_min = min(FQ);        % Minimum OFs

switch options.ranking

    case 'matlab' % --> do ranking via fast nondominated sorting in MATLAB
        Fr = {[]};              % Pareto-optimal fronts
        Sp = cell(M,1);         % Set points a particular point dominates
        np = zeros(M,1);        % # individuals which dominate a point
        % Now loop over all elements of R
        for p = 1 : M
            id = find((sum(bsxfun(@le,FQ(p,1:m),FQ),2)==m)& ...
                (sum(bsxfun(@lt,FQ(p,1:m),FQ),2)>0)); % Which points Q, p dominates
            if ~isempty(id), Sp{p} = id; end        % Store these in Sp
            id = find((sum(bsxfun(@le,FQ,FQ(p,1:m)),2)==m)& ...
                (sum(bsxfun(@lt,FQ,FQ(p,1:m)),2)>0)); % Which points Q dominate p
            np(p) = np(p) + numel(id);              % Store # points in np
            if np(p) == 0                           % p member 1st front
                Fr{1}(end+1) = p;
            end
        end
        i = 1;                                      % Determine rank
        while ~isempty(Fr{i})
            NextFr = [];                            % Next front is empty
            for p = 1 : numel(Fr{i})                % Each member p of Fr{i}
                q = Sp{Fr{i}(p)};                   % Modify each member Sp
                np(q) = np(q) - 1;                  % Decrease np(q) by one
                NextFr = [NextFr;q((np(q) == 0))];  %#ok if n(q) is zero, q is member of next front
            end
            i = i + 1;                              % Next front
            Fr{i} = NextFr;                         % Current front, members NextFr
        end
        for j = 1:i, RQ(Fr{j}) = j; end             % Extract N-vector with ranks from cell structure Fr

    case 'c' % --> now do ranking with help of C++
        for r = 1:M                                 % Loop unknown ranks points Q 
            id = find(isnan(RQ));                   % Which of Q no rank yet?
            if isempty(id)                          % break if idx empty 
                break 
            end
            front = paretofront(FQ(id,1:m));        % Find rank 1 solutions using C++ script
            ii = id(front == 1);                    % ID rank 1 idx batch?
            RQ(ii) = r;                             % Assign rank r 
        end

    otherwise
        error('do not know this language for ranking')
end

dQ = distance(FQ,RQ,options);

end

function dQ = distance(FQ,RQ,options)
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% %
% Fast nondominated sorting for Nxm matrix of objective function values   %
%                                                                         %
%  SYNOPSIS: rank = FNS(FX,RQ)                                            %
%   where                                                                 %
%    FQ        [input] Mxm matrix objective functions parents & children  %
%    RQ        [input] Mx1 vector of Pareto ranks N parents & children    %
%    options   [input] AMALGAM input structure algorithmic settings       %
%    .density       Density of rank 1,2,3, etc. solutions DEF: 'crowding' %
%     = 'crowding'  Crowding distance:Deb et al: NSGA-II  DEFault         %
%     = 'strength'  Strength Pareto:Zitzler&Thiele: SPEA-2                %
%    dQ        [outpt] Mx1 vector of crowding distance/strength Pareto's  %
%                                                                         %
%  © Written by Jasper A. Vrugt, Jan. 2005                                %
%  Los Alamos National Laboratory                                         %
%  University of California Irvine                                        %
%                                                                         %
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% %
[M,m] = size(FQ);
% Now determine the "distance" of each point of Q
switch options.density
    case 'crowding'                                 % Crowding distance Deb et al. (NSGA-II)
        dQ = nan(M,1);                              % Init. Mx1 vector crowding distances        
        for j = 1:max(RQ)
            id = (RQ == j);                         % Points with rank j
            R_s = FQ(id,1:m);                       % Select those points
            N_s = sum(id);                          % # points with rank j
            Cr_d = zeros(N_s,1);                    % (Re)-initialize crowding distance
            for i = 1:m
                [R_sort,sort_id] = sort(R_s(:,i));  % Sort OF ascending order
                Cr_d(sort_id([1 N_s])) = inf;       % Infinite distance "end" solutions
                                                    % Note: Cr_d(sort_idx(1),1) = Cr_d(sort_idx(1),1) + inf; 
                                                    %       is same as Cr_d(sort_idx(1),1) = inf;                      
                for z = 2:(N_s - 1)                 % Crowding distance others
                    Cr_d(sort_id(z)) = ...  
                       Cr_d(sort_id(z)) + ...
                       (R_sort(z+1) - R_sort(z-1));
                end
            end
            dQ(id,1) = Cr_d;                        % Crowding distance id solutions Q
        end
        
    case 'strength'                                 % Strength Pareto Zitlzer&Thiele (SPEA-2)
        n_dom = nan(M,1);                           % Init. M x 1 domination vector
        for j = 1 : M
            n_dom(j) = sum(sum(bsxfun(@lt, ...      % # points point Q dominates
                FQ(j,1:m),FQ),2)==m);       
        end
        dQ = M ./ n_dom;                            % Higher n_dom, lower strength (= many points nearby niche)
                                                    % Strength as proxy for crowding distance
    otherwise
        error('Do not know this density estimation method');
end

end
