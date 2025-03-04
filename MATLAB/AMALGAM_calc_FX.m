function [FX,Y] = AMALGAM_calc_FX(X,AMALGAMPar,f_handle,options,verbose)
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
% Evaluate user-supplied function and return objective function values    %
% model simulations (if so desired)                                       %
%                                                                         %
%  SYNOPSIS                                                               %
%   [FX,Y] = AMALGAM_calc_FX(X,AMALGAMPar,options,f_handle,verbose)       %
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

if nargin < 5
    verbose = 0;
else
    clear('prt_progress');      % Waitbar initialization
end
[N,d] = size(X);                % # parameter vectors?
X = X';                         % Transpose matrix X (column vectors)
m = AMALGAMPar.m;               % # objective functions
FX = nan(N,m); Y = [];          % Preallocate FX and model simulations

% Evaluate model
switch AMALGAMPar.CPU                                 
    case 1    % Sequential evaluation
        switch options.modout
            case 'no'
                for ii = 1:N
                    FX(ii,1:m) = f_handle(X(1:d,ii)); 
                    if verbose, prt_progress(AMALGAMPar,N); end
                end
            case 'yes'
                % Execute the model and return the objective functions
                for ii = 1:N
                    [FX(ii,1:m),Y(ii,:)] = f_handle(X(1:d,ii));     %#ok
                    if verbose, prt_progress(AMALGAMPar,N); end
                end
        end
    otherwise   % Parallel evaluation
        if strcmp(options.IO,'yes')
            EXAMPLE_dir = pwd;
            % Model call depends on whether we return model simulation or not
            switch options.modout
                case 'no'
                    parfor ii = 1:N
                        % Determine work ID
                        task = getCurrentTask(); id = get(task, 'ID');
                        % Go to right directory (t.id is directory number)
                        cd(strcat(EXAMPLE_dir,slash_dir,num2str(id)));
                        % Execute model, return objective functions
                        FX(ii,:) = f_handle(X(1:d,ii));             %#ok
                        % Print output
                        if verbose && id == 1
                            prt_progress(AMALGAMPar,N);
                        end
                    end
                case 'yes'
                    parfor ii = 1:N
                        % Determine work ID
                        task = getCurrentTask(); id = get(task, 'ID');
                        % Go to right directory (t.id is directory number)
                        cd([EXAMPLE_dir,'/',num2str(id)]);
                        % Execute model, objective functions & simulations
                        [FX(ii,:),Y(ii,:)] = f_handle(X(1:d,ii));   %#ok
                        % Print output
                        if verbose && id == 1
                            prt_progress(AMALGAMPar,N);
                        end
                    end
            end
            cd(EXAMPLE_dir)
        elseif strcmp(options.IO,'no')
            switch options.modout
                case 'no'
                    parfor ii = 1:N
                        % Determine work ID
                        task = getCurrentTask(); id = get(task, 'ID');
                        % Execute model, return objective functions
                        FX(ii,:) = f_handle(X(1:d,ii));             %#ok
                        % Print output
                        if verbose && id == 1
                            prt_progress(AMALGAMPar,N);
                        end
                    end
                case 'yes'
                    parfor ii = 1:N
                        % Determine work ID
                        task = getCurrentTask(); id = get(task, 'ID');
                        % Execute model, objective functions & simulations
                        [FX(ii,:),Y(ii,:)] = f_handle(X(1:d,ii));   %#ok
                        % Print output
                        if verbose && id == 1
                            prt_progress(AMALGAMPar,N);
                        end
                        
                    end
            end
        end
end

% Go to next line if verbose is activated
if verbose, fprintf('\n'); fprintf('Model simulation ... done\n'); end

end
