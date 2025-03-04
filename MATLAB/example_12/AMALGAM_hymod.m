function [F,Y_sim] = AMALGAM_hymod(par,plugin)
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% %
% hymod simulation of discharge and returns a pair of objective function  %
% values, namely, the Root Mean Square Errors of driven and nondriven     %
% parts of the hydrograph                                                 %
%                                                                         %
%  REFERENCES                                                             %
%  Vrugt, J.A., H.V. Gupta, L.A. Bastidas, W. Bouten, and S. Sorooshian   %
%      (2000), Effective and efficient algorithm for multiobjective       %
%      optimization of hydrologic models, Water Resources Research,       %
%      39 (8), 1214, https://doi.org./10.1029/2002WR001746                %
%                                                                         %
%  © Written by Jasper A. Vrugt, Jan. 2006                                %
%  Los Alamos National Laboratory                                         %
%  University of California Irvine                                        %
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% %

p = num2cell(par); [cmax,bexp,alfa,Rs,Rq] = ... % Extract parmeters
    deal(p{:}); 
[T_max,Y_obs,PET,R,idx_d,N_d,idx_nd,N_nd] = ... % Extract data
    v2struct(plugin,plugin.fields);

x_l = 0.0;                                      % Initial state loss tank
x_s = 0;                                        % Initial state slow tank
x_q = zeros(3,1);                               % Initial states fast tanks
output = nan(T_max,1);                          % Simulated discharge

for t = 1 : T_max   % Now loop over the forcing data
    % Compute excess precipitation and evaporation
    [ER1,ER2,x_l] = excess(x_l,cmax,bexp,R(t,1),PET(t,1));
    % Calculate total effective rainfall
    ET = ER1 + ER2;
    % Now partition ER between quick and slow flow reservoirs
    UQ = alfa * ET; US = (1-alfa) * ET;
    % Route slow flow component with single linear reservoir
    [x_s,QS] = linres(x_s,US,Rs);
    % Route quick flow component with linear reservoirs
    q_in = UQ;
    % Series of three linear fast reservoirs
    for k = 1 : 3
        [x_q(k),q_out] = linres(x_q(k),q_in,Rq); 
        q_in = q_out;
    end
    % Compute total flow for timestep (in mm/day)
    output(t,1) = (QS + q_out); 
end
Y_sim = output(65:T_max,1);                     % Apply burn-in of 65 days
F(1,1) = sqrt(sum((Y_sim(idx_d) - ...           % RMSE driven part 
    Y_obs(idx_d)).^2)/N_d);
F(2,1) = sqrt(sum((Y_sim(idx_nd) - ...          % RMSE nondriven part 
    Y_obs(idx_nd)).^2)/N_nd);

end
%<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>

%<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>
% Secondary functions
%<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>

% 1: NSGA_crossover
function [ER1,ER2,xn] = excess(x_loss,cmax,bexp,Pval,PETval)
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% %
%excess: Function calculates excess precipitation and evaporation
%
% Written by Jasper A. Vrugt
% Los Alamos National Laboratory
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% %

xn_prev = x_loss;
ct_prev = cmax*(1-power((1-((bexp+1)*(xn_prev)/cmax)),(1/(bexp+1))));
% Calculate effective rainfall 1
ER1 = max((Pval-cmax+ct_prev),0.0);
Pval = Pval-ER1;
dummy = min(((ct_prev+Pval)/cmax),1);
xn = (cmax/(bexp+1))*(1-power((1-dummy),(bexp+1)));
% Calculate effective rainfall 2
ER2 = max(Pval-(xn-xn_prev),0);

% Alternative approach
evap = (1-(((cmax/(bexp+1))-xn)/(cmax/(bexp+1))))*PETval; 
% Actual ET is linearly related to the soil moisture state
xn = max(xn-evap, 0);

%evap = min(xn,PETval);
%xn = xn-evap;

end
%<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>

% 2: linres
function [x,q_out] = linres(x,q_in,R)
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% %
%linres: Function of linear reservoir
%
% Written by Jasper A. Vrugt
% Los Alamos National Laboratory
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% %

x = (1-R)*x + (1-R)*q_in;
q_out = (R/(1-R))*x;

end
%<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>
