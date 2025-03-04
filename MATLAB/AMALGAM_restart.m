function [AMALGAMPar,Par_info,options,PS,f_handle,X,Z,FX,YX,p_rm,output, ...
    FX_min,RX,dX,ct,T_start] = AMALGAM_restart(fname,Ftrue)            %#ok
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
% Restart function to complete desired number of generations              %
%                                                                         %
%  SYNOPSIS                                                               %
%   [AMALGAMPar,Par_info,options,PS,f_handle,X,Z,FX,YX,p_rm,output, ...   %
%       FX_min,RX,dX,ct,T_start] = AMALGAM_restart(fname,Ftrue)           %
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

% Try ... catch statement in case restart .mat file does not exist
try
    % String with data to read from fname
    evalstr = strcat('load',{' '},fname,{' '},['AMALGAMPar Par_info ' ...
        ['options PS X Z FX YX output FX_min RX dX p_rm' ...
        'Func_name plugin Func_name t']]);
    % Load the restart data
    eval(char(evalstr));
catch
    % Create error warning to screen
    evalstr = char(strcat(['DREAM_PACKAGE ERROR: Cannot restart --> ' ...
        'File '],{' '},fname,{' '},['does not exist. Next run, to ' ...
        'avoid this problem, set field ''save'' of structure options' ...
        ' to ''yes'' ']));
    % Now print error
    error(evalstr);
end

% Open warning file
fid = fopen('warning_file.txt','a+');

% Check whether previous run was aborted early or not
switch ( t < DREAMPar.T )
    case 1      % Finish previous budget
        evalstr = char(strcat(['AMALGAM RESTART: Starting with ' ...
            't = '],{' '},num2str(t),{' '},['but still using old ' ...
            'budget of'],{' '},'AMALGAMPar.T = ',{' '}, ...
            num2str(AMALGAMPar.T),'\n'));                           %#ok
        % Now print warning to screen and to file
        fprintf(evalstr); fprintf(fid,evalstr);

    otherwise   % --> assign new budget and run
        fprintf(fid,'-------------- AMALGAM warning file --------------\n');
        fname = 'T.txt';    % Now load file
        % Check whether file T.txt exists
        if exist(fname,'file') == 2 % Yes
            evalstr = char(strcat(['AMALGAM RESTART: Located ' ...
                'file ''T.txt'' in respective example ' ...
                'directory - now checking its content\n']));
            % Now print warning to screen and to file
            fprintf(evalstr); fprintf(fid,evalstr);
            % Now check content of T.txt to avoid crash of code
            T_new = checkfile_T_AMALGAM(fname);
            % If all are satisfied then
            evalstr = char(strcat(['AMALGAM RESTART: User has ' ...
                'requested/listed'],{' '},num2str(T_new),{' '}, ...
                'additional generations in file ''T.txt''\n'));
            % Now print warning to screen and to file
            fprintf(evalstr); fprintf(fid,evalstr);
            % Now print to file/screen the final message
            evalstr = char(strcat('AMALGAM RESTART: Initial t = ', ...
                {' '},num2str(t),{' '},'and completing',{' '}, ...
                num2str(T_new),{' '},...
                'additional generations so that',{' '},'AMALGAMPar.T = ', ...
                {' '},num2str(AMALGAMPar.T  + T_new),'\n'));        %#ok
            % Now print warning to screen and to file
            fprintf(evalstr); fprintf(fid,evalstr);

        else % file does not exist
            evalstr = char(strcat(['AMALGAM RESTART: Did not locate ' ...
                'file ''T.txt'' in respective example directory so ' ...
                'per default we double the "previous" budget of ' ...
                'generations\n']));
            % Now print warning to screen and to file
            fprintf(evalstr); fprintf(fid,evalstr);
            % Thus
            T_new = AMALGAMPar.T;                                   %#ok
            % Now print to file/screen the final message
            evalstr = char(strcat('AMALGAM RESTART: Initial t = ',{' '}, ...
                num2str(t),{' '},'and completing',{' '},num2str(T_new), ...
                {' '},'additional generations so that',{' '}, ...
                'AMALGAMPar.T = ',{' '}, ...
                num2str(AMALGAMPar.T  + T_new),'\n'));
            % Now print warning to screen and to file
            fprintf(evalstr); fprintf(fid,evalstr);
        end
        [T,d2] = size(Z);                                           %#ok                          
        % Add nan to Z
        Z(T+1:T+T_new,1:d2) = nan(T_new,d2);
        % Update DREAMPar.T with T_new
        AMALGAMPar.T = AMALGAMPar.T + T_new;
end

% close warning_file.txt
fclose(fid);
% Define starting value of T
T_start = t + 1;
% Setup parallel computing framework or not
[AMALGAMPar,f_handle] = AMALGAM_calc_setup(AMALGAMPar,Func_name, ...
    options,plugin);

end