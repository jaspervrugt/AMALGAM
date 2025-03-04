function [AMALGAMPar,Par_info,options,T_start] = ...
    AMALGAM_setup(AMALGAMPar,Par_info,options)
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
% Initializes the main variables used in AMALGAM                          %
%                                                                         %
%  SYNOPSIS                                                               %
%   [AMALGAMPar,Par_info,options,T_start] = AMALGAM_setup( ...            %
%       AMALGAMPar,Par_info,options)                                      %
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

% Random seed (legacy: randn('state',sum(100*clock)); )
rng(1+round(100*rand),'twister');

% Define T_start
T_start = 2;

% Field names of AMALGAMPar
f_names = {'beta_1','beta_2','c_1','c_2','varphi','p_CR', ...
    'p_M','eta_C','eta_M','gamma','K','p0'};
% Default settings
value = {'@(N) unifrnd(0.6,1.0,N,1)','@(N) unifrnd(0.2,0.6,N,1)', ...
    '1.5','1.5','@(N) unifrnd(0.5,1.0,N,1)','0.9','1/AMALGAMPar.d', ...
    '20','20','(2.38/sqrt(AMALGAMPar.d))^2','1','0.05'};

% Now check
for j = 1:numel(f_names)
    if ~isfield(AMALGAMPar,f_names(j))
        AMALGAMPar.(char(f_names(j))) = eval(char(value(j)));
    end
end

% If recombination methods not specified --> use default AMALGAM
if ~isfield(AMALGAMPar,'rec_methods')
    AMALGAMPar.rec_methods = {'ga','ps','am','de'};
end

% Field names of options
f_names = {'parallel','IO','save','restart','modout','ranking', ...
    'density','print'};
% Default values/settings
value = {'no','no','no','no','no','matlab','crowding','yes'};
% Undefined fields, set to default
for j = 1 : numel(f_names)
    if ~isfield(options,f_names(j))
        options.(char(f_names(j))) = char(value(j));
    end
end

% Replicate Par_info.min
% Par_info.min = repmat(Par_info.min,AMALGAMPar.N,1); 
% Replicate Par_info.max
% Par_info.max = repmat(Par_info.max,AMALGAMPar.N,1); 
% Compute step size and make N copies
if isfield(Par_info,'steps')
    % Par_info.steps = repmat(Par_info.steps,AMALGAMPar.N,1);
    Par_info.step_size = (Par_info.max - Par_info.min) ./ Par_info.steps;
end

% Now check for BMA
global BMA;                                                         %#ok
if ~isempty(BMA)
    Par_info.unit_simplex = 'yes'; 
end

% Define prior handle (random samples & evaluation of pdf)
if isfield(Par_info,'prior')
    if iscell(Par_info.prior)
        Par_info.u = 'yes';                 % Univariate case
        for ii = 1 : AMALGAMPar.d           % Anonymous handle each paramtr
            Par_info.prior_rnd{ii} = ...    % Handle draw initial state
                eval(char(strcat('@(x)',{' '}, ...
                char(strrep(Par_info.prior(ii), ...
                'pdf','rnd')))));
        end
    else
        Par_info.u = 'no';                  % Multivariate case
        pr_name = char(Par_info.prior);     % Turn handle to string
        [pr_var,idcp] = ...                 % Extract variable names prior
            extract_names(pr_name);         % & index closing parenthesis
        pr_name = pr_name(idcp(1)+1:end);   % Prior handle without @()
        n_var = numel(pr_var);              % # variabls without x (=n_p-1)
                                            % n_p = nargin(Par_info.prior)
        for z = 1:n_var                     % Add Par_info. to varble names        
            pr_name = strrep(pr_name,...
                pr_var{z},strcat('Par_info.',pr_var{z}));
        end
        Par_info.prior_rnd = ...            % Handle draw initial states
            eval(char(strcat('@(x)',{' '}, ...
                char(strrep(pr_name, ...
                'pdf(x,','rnd(')))));
    end
end

% Now print to screen AMALGAMPar fields
AMALGAMPar = orderfields(AMALGAMPar,{'d','N','T','m','rec_methods', ...
    'K','p0','beta_1','beta_2','c_1','c_2','varphi','p_CR','p_M', ...
    'eta_C','eta_M','gamma'});
not_print = 0;

% Now print to screen
fprintf('--- Summary of algorithmic settings ---\n');
fprintf('  AMALGAMPar\n');
f_names = fieldnames(AMALGAMPar);
for i = 1:numel(f_names) - not_print
    pr = eval(char(strcat('AMALGAMPar.',f_names(i,:))));
    if isreal(pr)
        fprintf(char(strcat('    .',f_names(i,:),{' '},'=',{' '}, ...
            num2str(pr),'\n')));
    elseif iscell(pr)
        pr = char(pr); prt = strcat('{',pr(1,:)); 
        for u = 2:size(pr,1)
            prt = strcat(prt,',',pr(u,:)); 
        end 
        prt = strcat(prt,'}');
        fprintf(char(strcat('    .',f_names(i,:),{' '},'=',{' '}, ...
            prt,'\n')));
    else
        fprintf(char(strcat('    .',f_names(i,:),{' '},'=',{' '}, ...
            func2str(pr),'\n')));
    end
end

fprintf('----- End of algorithmic settings -----\n');

end
