function [AMALGAMPar,f_handle] = AMALGAM_calc_setup(AMALGAMPar, ...
    fname,options,plugin)
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
% Sets up sequential / parallel                                           %
%                                                                         %
%  SYNOPSIS                                                               %
%   [AMALGAMPar,f_handle] = AMALGAM_calc_setup(AMALGAMPar,Func_name, ...  %
%       options,plugin)                                                   %
%  where                                                                  %
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

if isempty(plugin)
    f_handle = eval(['@(x)',char(fname),'(x)']);        % Function handle
else
    f_handle = eval(['@(x)',char(fname),'(x,plugin)']);
end
T = 2;                                                  % First generation
switch options.parallel                                 % Parallel chain execution?
    case 'no'
        AMALGAMPar.CPU = 1;                             % We use 1 CPU (processor)
    case 'yes'
        if verLessThan('matlab','8.2')                  % MATLAB version
            vMATLAB = 'old';
        else
            vMATLAB = 'new';
        end       
        switch vMATLAB                                  % Open cluster, MATLAB version
            case 'old'  % Use matlabpool
                isOpen = matlabpool('size') > 0;        %#ok<DPOOL> Cluster open?
                if ~isOpen, matlabpool open, end        %#ok<DPOOL> If not, new pool created 
                workers = matlabpool('size');           %#ok<DPOOL> # processors
                if ( workers > AMALGAMPar.N )           % Close the cluster
                    matlabpool close force local        %#ok<DPOOL> 
                    evalstr = strcat(['matlabpool ' ... % Reopen DREAMPar.N cores
                        'open'],{' '}, ...
                        num2str(AMALGAMPar.N)'); 
                    eval(char(evalstr));
                    workers = AMALGAMPar.N;             % # workers is DREAMPar.N
                end
            case 'new'  % Use parpool
                pool = gcp;                             % Cluster open?
                if isempty(pool)                        % If not, new pool created
                    pool = parpool('local'); 
                end
                workers = pool.NumWorkers;              % # processors
                if ( workers > AMALGAMPar.N )
                    delete(gcp)                         % Close the cluster
                    parpool('local',AMALGAMPar.N);      % Reopen cluster DREAMPar.N cores
                    workers = AMALGAMPar.N;             % # workers is DREAMPar.N
                end
        end
        AMALGAMPar.CPU = workers;                       % DREAMPar.CPU is "workers"
        fid = fopen('warning_file.txt','a+');           % Open warning_file.txt
        evalstr = char(strcat(['AMALGAM ' ...           % Print to screen
            'PARALLEL: MATLAB pool ' ...                
            'opened with'],{' '},num2str(AMALGAMPar.CPU), ...
            {' '},'workers for',{' '}, ...
            num2str(AMALGAMPar.N),{' '},'chains \n'));
        fprintf(evalstr); fprintf(fid,evalstr);         % Print warning to screen & file
        fclose(fid);                                    % Close warning.txt
    
        if strcmp(options.IO,'yes')                     % If input/output writing, 
                                                        % we need directories each worker
            if (ispc || ismac)                          % Copy files to work directories
                copyfile('*.*',[pwd,'/temp_AMALGAM'])   % copy first to temporary directory
                for ii = 1:AMALGAMPar.CPU               % now move to this directory
                    copyfile([pwd,'/temp_AMALGAM/', ...
                        '*.*'],[pwd,'/',num2str(ii)]);
                end
                rmdir([pwd,'/temp_AMALGAM'],'s');       % remove temporary directory
            end
            % LINUX/UNIX ENVIRONMENT BUT WORKS FOR PC AS WELL
        end
    
end
    
end
