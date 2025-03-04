function AMALGAM_postproc_late
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
% Function postprocesses output of AMALGAM outside main PACKAGE           %
%                                                                         %
%  SYNOPSIS                                                               %
%   AMALGAM_postproc_late                                                 %
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

fname = 'AMALGAM.mat';                          % File name of memory
load_str = strcat('load',{' '},fname, ...       % Create load string    
    {' '},['AMALGAMPar Par_info ' ...
    ['options X FX output Ftrue YX Z ct t0']]);
eval(char(load_str));                           % Load the data

output.RunTime = cputime - t0;                  % Determine total run time
AMALGAM_end(AMALGAMPar,options);                % Close workers
Z = Z(1:ct*AMALGAMPar.N,1:AMALGAMPar.d ...      % Finalize achive Z
    + AMALGAMPar.m);
AMALGAM_postproc(AMALGAMPar,Par_info, ...       % Postprocess results
    options,X,FX,output,Ftrue,YX,Z);

% Return arguments: X FX output Z YX

end