function [F,Y_sim] = AMALGAM_hmodel(par,plugin)
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% %
% hmodel simulation of discharge according to Schoups et al. 2010         %
% Returns Root Mean Square Error of driven and nondriven part hydrograph  %
%                                                                         %
%  REFERENCES                                                             %
%  Schoups, G., and J.A. Vrugt (2010), A formal likelihood function for   %
%      parameter and predictive inference of hydrologic models with       %
%      correlated, heteroscedastic, and non-Gaussian errors, Water        %
%      Resources Research, 46, W10531,                                    %
%      https://doi.org/10.1029/2009WR008933                               %
%  Vrugt, J.A., H.V. Gupta, L.A. Bastidas, W. Bouten, and S. Sorooshian   %
%      (2000), Effective and efficient algorithm for multiobjective       %
%      optimization of hydrologic models, Water Resources Research,       %
%      39 (8), 1214, https://doi.org./10.1029/2002WR001746                %
%                                                                         %
%  © Written by Jasper A. Vrugt, Jan. 2011                                %
%  Los Alamos National Laboratory                                         %
%  University of California Irvine                                        %
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% %
 
[tout,data,hmodel_opt,y0,Y_obs,n,id_d, ...  % Extract various input data
    N_d,id_nd,N_nd] = ...
    v2struct(plugin,plugin.fields);
y = hmodel(par,tout,data,hmodel_opt,y0);    % Simulate discharge with par
Y = y(5,2:n+1) - y(5,1:n);                  % Compute discharge from state
Y_sim = Y(731:n);                           % Apply burn-in [= 2 years]
F(1) = sqrt(sum((Y_sim(id_d) - ...          % RMSE driven part hydrograph
    Y_obs(id_d)).^2)/N_d);     
F(2) = sqrt(sum((Y_sim(id_nd) - ...         % RMSE nondrivn part hydrograph
    Y_obs(id_nd)).^2)/N_nd);  

end
