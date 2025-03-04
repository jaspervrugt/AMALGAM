function AMALGAM_end(AMALGAMPar,options)
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
% Now close workers and delete worker directories, if appropriate         %
%                                                                         %
%  SYNOPSIS                                                               %
%   AMALGAM_end(AMALGAMPar,options)                                       %
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

if AMALGAMPar.CPU > 1
    if verLessThan('matlab','8.2')      
        parpool('close');               % Close MATLAB pool
    else
        delete(gcp('nocreate'));        % Close in latest version
    end
    if strcmp(options.IO,'yes')         % If IO writing, remove directories
        try                             % Remove each worker directory
            for ii = 1:AMALGAMPar.CPU
                rmdir([pwd,'/',num2str(ii)],'s');
            end
        catch                           % Otherwise, return to header
            return 
        end
    end
end

% Open warning_file.txt
fid = fopen('warning_file.txt','a+','n');

% Write final line of warning file
fprintf(fid,'----------- End of AMALGAM warning file ----------\n');
% Close the warning file
fclose(fid);

% Open the warning file
if ispc || ismac, edit warning_file.txt, end

% % % Now create a string for plotting of results
% % str = strcat('{','''x_{1}'',');
% % for i = 2 : AMALGAMPar.d - 1
% %     str = strcat(str,'''x_{',num2str(i),'}'',');
% % end
% % 
% % % Create string for printing of weights to screen
% % str = strcat(str,'''x_{',num2str(AMALGAMPar.d),'}''}'); str = eval(str);
% % 
% % % Now create another string for plotting
% % str_plot = str; str_plot = strrep(str_plot,'x','$x'); 
% % str_plot = strrep(str_plot,'}','}$');   
% % 
% % % Now store for plotting tables and figures
% % plotting.str_plot = str_plot; plotting.str = str; 