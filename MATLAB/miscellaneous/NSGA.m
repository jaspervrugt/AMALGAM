function G = NSGA(AMALGAMPar,Par_info,X,FX,id,dX)
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% %
% Generate offspring using NSGA-II method                                 %
%                                                                         %
%  SYNOPSIS                                                               %
%   G = NSGA(AMALGAMPar,Par_info,X,F,id,dX)                               %
%  where                                                                  %
%                                                                         %
%  © Written by Jasper A. Vrugt, Jan. 2005                                %
%  Los Alamos National Laboratory                                         %
%  University of California Irvine                                        %
%                                                                         %
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% %

G = NSGA_crossover(AMALGAMPar,Par_info,X,FX,dX);% Selection & crossover
G = NSGA_mutate(AMALGAMPar,Par_info,G);         % Polynomial mutation
G = G(id,1:AMALGAMPar.d);                       % Select individuals

end
%<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>

%<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>
% Secondary functions
%<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>

% 1: NSGA_crossover
function G = NSGA_crossover(AMALGAMPar,Par_info,X,FX,dX)
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% %
%NSGA_crossover: Function performs genetic selection and crossover
%
% Written by Jasper A. Vrugt based on C++ code Deb
% Los Alamos National Laboratory
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% %

% If population size is uneven - need to account for this with offspring
div = AMALGAMPar.N/2; 
% Check whether population size is even or uneven
res = ceil ( div - floor(div) );
% Define a1 and a2
a1 = randperm(AMALGAMPar.N); a2 = randperm(AMALGAMPar.N);
% Adjust a1 and a2 if uneven (then length becomes even)
a1 = [a1 a1(1:res)]; a2 = [a2 a2(1:res)]; 
% Initialize the new Children Generation
G = []; ct = 1;
% Now do tournament selection and crossover
while size(G,1) < AMALGAMPar.N   
    % Select parents
    for j = 1 : 2
        if j == 1
            a_1 = a1(ct); a_2 = a1(ct + 1); 
            a_3 = a1(ct + 2); a_4 = a1(ct + 3);
        else
            a_1 = a2(ct); a_2 = a2(ct + 1); 
            a_3 = a2(ct + 2); a_4 = a2(ct + 3);
        end
        % Tournament selection
        parent1 = Tournament(AMALGAMPar,X(a_1,:),X(a_2,:), ...
            FX(a_1,:),FX(a_2,:),dX(a_1),dX(a_2));
        % Tournament selection
        parent2 = Tournament(AMALGAMPar,X(a_3,:),X(a_4,:), ...
            FX(a_3,:),FX(a_4,:),dX(a_3),dX(a_4));
        % Now do crossover
        [child1,child2] = crossover(AMALGAMPar,Par_info,parent1,parent2);
        % Append the children to the offspring
        G = [G ; child1 ; child2];                                  %#ok
        % Update counter
        ct = ct + 1;
    end
end
% Select the first AMALGAMPar.N individuals
G = G(1:AMALGAMPar.N,1:AMALGAMPar.d);

end
%<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>

% 2: Tournament
function [parent] = Tournament(AMALGAMPar,a,b,Fa,Fb,da,db)
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% %
%Tournament: Function performs tournament selection
%
% Written by Jasper A. Vrugt based on C++ code Deb
% Los Alamos National Laboratory
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% %

% First draw random number
rnd = rand; done = 0;
% First check dominance
flagout = Check_dominance(AMALGAMPar,Fa,Fb);
% Now check various flags
if flagout == 1
    parent = a; done = 1;
end
if flagout == -1
    parent = b; done = 1;
end
if (da > db) && (done == 0)
    parent = a; done = 1;
end
if (da < db) && (done == 0)
    parent = b; done = 1;
end
if (rnd <= 0.5) && (done == 0)
    parent = a;
end
if (rnd > 0.5) && (done == 0)
    parent = b;
end

end
%<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>

% 3: Check_dominance
function [flagout] = Check_dominance(AMALGAMPar,Fa,Fb)
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% %
%Tournament: Function checks for dominance [multiobjective sense]
%
% Written by Jasper A. Vrugt based on C++ code Deb
% Los Alamos National Laboratory
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% %

% First initialize some important variables
flag1 = 0; flag2 = 0;
% Check objective function values
for qq = 1:AMALGAMPar.m
    if Fa(qq) < Fb(qq)
        flag1 = 1;
    end
    if Fa(qq) > Fb(qq)
        flag2 = 1;
    end
end
% Calculate flagout
delta = flag1 - flag2; flagout = -1;

switch delta
    case 0
        flagout = 0;
    case 1
        flagout = 1;
end

end
%<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>

% 4: crossover
function [g1,g2] = crossover(AMALGAMPar,Par_info,x1,x2)
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% %
%crossover: Function performs crossover of two individuals
%
% Written by Jasper A. Vrugt based on C++ code Deb
% Los Alamos National Laboratory
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% %

if rand < AMALGAMPar.p_CR
    % Define lower and upper bound value
    yl = Par_info.min(1,1:AMALGAMPar.d); 
    yu = Par_info.max(1,1:AMALGAMPar.d);
    % Create random numbers
    rnd = rand(1,AMALGAMPar.d); idx = find(rnd > 0.5);
    % Determine idx_par
    idx_par = find( abs(x1 - x2) <= eps);
    % Derive y1 and y2
    y1 = min(x1,x2); y2 = max(x1,x2);
    % Determine beta and alpha
    beta = 1.0 + (2.0*(y1-yl)./(y2-y1)); alpha = 2.0 - ...
        beta.^(-(AMALGAMPar.eta_C+1.0));
    % Create random number
    rnd = rand(1,AMALGAMPar.d);
    % Determine ii1 and ii2
    ii1 = find(rnd <= (1./alpha)); ii2 = find(rnd > (1./alpha));
    % Prior define betaq
    betaq = 1e-10 * ones(1,AMALGAMPar.d);
    % Determine betaq (ii1)
    betaq(1,ii1) = (rnd(ii1).*alpha(ii1)).^(1.0/(AMALGAMPar.eta_C+1.0));
    % Determine betaq (ii2)
    betaq(1,ii2) = (1.0./(2.0 - rnd(ii2).*alpha(ii2))).^ ...
        (1.0/(AMALGAMPar.eta_C+1.0));
    % Determine g1
    g1 = 0.5 * ( ( y1 + y2 ) - betaq.* (y2 - y1 ) );
    % And determine beta and alpha again
    beta = 1.0 + ( 2.0 * ( yu - y2 )./ ( y2 - y1 ) ); 
    alpha = 2.0 - beta.^(-(AMALGAMPar.eta_C+1.0));
    % Determine betaq (ii1)
    betaq(1,ii1) = (rnd(ii1).*alpha(ii1)).^ ...
        ( 1.0 / (AMALGAMPar.eta_C + 1.0 ) );
    % Determine betaq (ii2)
    betaq(1,ii2) = (1.0./(2.0 - rnd(ii2).*alpha(ii2))).^ ...
        (1.0/(AMALGAMPar.eta_C+1.0));
    % Determine g2
    g2 = 0.5*((y1+y2)+betaq.*(y2-y1));
    % Make sure that g1 is within bound
    g1 = max(g1,yl); g1 = min(g1,yu); g2 = max(g2,yl); g2 = min(g2,yu);
    % Create random number
    rnd = rand(1,AMALGAMPar.d);
    % Now determine which values of c1 to use
    ii1 = find(rnd <= 1/2); ii2 = find(rnd > 1/2);
    % And create child 1 and 2
    g1(1,ii1) = g2(ii1); % g2(1,ii1) = g1(ii1);
    % And create child 1 and 2
    g1(1,ii2) = g1(ii2); % g2(1,ii2) = g2(ii2);
    % change if ( abs ( parent1 - parent2 ) < eps )
    g1(1,idx_par) = x1(1,idx_par); g2(1,idx_par) = x2(1,idx_par);
    % Now update again (50% chance to replace with parent)
    g1(1,idx) = x1(1,idx); g2(1,idx) = x2(1,idx);
else    % No crossover
    g1 = x1; g2 = x2;
end

end
%<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>

% 5: NSGA_mutate
function G = NSGA_mutate(AMALGAMPar,Par_info,G)
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% %
%NSGA_mutate: Function performs polynomial mutation
%
% Written by Jasper A. Vrugt based on C++ code Deb
% Los Alamos National Laboratory
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% %

% Calculate delta1
delta1 = (G - Par_info.min)./( Par_info.max - Par_info.min );
% Calculate delta2
delta2 = (Par_info.max - G)./( Par_info.max - Par_info.min );
% Draw uniform random number
rnd = rand ( AMALGAMPar.N , AMALGAMPar.d );
% Compute mutation power
mut_pow = 1/(AMALGAMPar.eta_M+1);
% Which values leq 0.5
idx = find(rnd <= 1/2); 
% Calculate xy and val
xy(idx,1) = 1 - delta1(idx); val(idx,1) = 2 * rnd(idx) + ...
    (1 - 2 * rnd(idx)).*(xy(idx,1).^(AMALGAMPar.eta_M + 1));
% Calculate deltaq
deltaq(idx) = val(idx,1).^mut_pow - 1;
% Which values g 0.5
idx = find(rnd > 1/2);
% Calculate xy and val
xy(idx,1) = 1 - delta2(idx); val(idx,1) = 2 * ( 1 - rnd(idx) ) + ...
    2 * (rnd(idx) - 1/2).*(xy(idx,1).^(AMALGAMPar.eta_M + 1));
% Calculate deltaq
deltaq(idx) = 1 - val(idx,1).^mut_pow;
% Check found bounds
z = G + reshape ( deltaq , AMALGAMPar.N , AMALGAMPar.d ) .* ...
    ( Par_info.max - Par_info.min );
% Make sure proposals are within bound
% z = max(z,Par_info.min); z = min(z,Par_info.max);
% Draw uniform random number
rnd = rand ( AMALGAMPar.N , AMALGAMPar.d );
% Now determine which values of Y are mutated
idx = find ( rnd <= AMALGAMPar.p_M );
% Mute those coordinates of Y
G(idx) = z(idx);

end
%<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>

