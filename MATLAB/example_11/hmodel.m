function [y] = hmodel(x,tout,data,options,y0)
% Runs the hmodel and returns the driven and nondriven part
%
% This code runs a "C" code using mex compiler in MATLAB. Please compile
% the C code using the mex compiler. You only need to do this one time.
% If the code is not working please contact me directly. The C-code and mex
% file is developed in MATLAB 2010, and might create problems with newer
% MATLAB versions.

%% Generate mex file (need to do this only once)
% mex crr_model.c;

data.Imax  = x(1);          % interception storage capacity (mm)
data.Sumax = x(2);          % unsaturated zone storage capacity (mm)
data.Qsmax = x(3);          % maximum percolation rate (mm/d)
data.aE    = x(4);          % evaporation coefficient
data.aF    = x(5);          % runoff coefficient
data.aS    = 1e-6;          % percolation coefficient
data.Kf    = x(6);          % fast-flow response time (d)
data.Ks    = x(7);          % slow-flow response time (d)

%% Run model C
y = crr_model(tout,y0,data,options);