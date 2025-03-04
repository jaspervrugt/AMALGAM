function AMALGAM_postproc(AMALGAMPar,Par_info,options,X,FX, ...
    output,Fpareto,YX,Z)
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
% Function postprocesses output of AMALGAM and creates tables & figures   %
%                                                                         %
%  SYNOPSIS                                                               %
%   AMALGAM_postproc(AMALGAMPar,Par_info,options,X,FX,output, ...         %
%       Fpareto,YX,Z)                                                     %
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

% Set sim to be empty if not specified
if nargin < 8
    YX = [];
end

set(0,'defaultTextInterpreter','latex'); %trying to set the default
colororder =  [ 0   0   1       % 1 BLUE
    1   0   1       % 5 MAGENTA (pale)
    0   1   0       % 2 GREEN (pale)
    1   .5  .25     % 11 ORANGE
    0   1   1       % 4 CYAN
    .6  .5  .4      % 15 BROWN (dark)
    0   .5  0       % 9 GREEN (dark)
    0   0   0 ];    % 7 BLACK
%n_colors = size(colororder,1);

% import mlreportgen.dom.*;
% doc = Document('HDMR_figures','pdf'); open(doc)
% Determine screen size
scr_z = get(0,'ScreenSize');
% Multiplier, x and y: axis
x_mult = scr_z(3)/1920; y_mult = scr_z(4)/1080;
% Multiplier, text
t_mult = min(x_mult,y_mult);
% Define fontsize for figures
fontsize_title = 20 * t_mult;
fontsize_axis_numbers = 20 * t_mult;
fontsize_labels = 20 * t_mult;
fontsize_A = 18 * t_mult;
linewidth_marker = 2 * t_mult;
fontsize_titlepage = 30 * t_mult;
fontsize_xylabel = 20 * t_mult;
fontsize_legend = 20 * t_mult;
% Define colors for different single criterion optima
symbol = {'r','b','g','c','k','m','y'};
% Associated colors
symbol_text = {'red','blue','green','cyan','black','magenta','yellow'};
% Maximum number of bins for histograms
maxbins = 25;

% Define name of program
n_program = 'AMALGAM';
% Define name of figures file
file_name = [n_program,'_figures.pdf'];

% Create legend/label string for different parameters
str_par = cell(AMALGAMPar.d,1);
switch isfield(Par_info,'names')
    case 1 % User specified the names of the parameters
        for i = 1:AMALGAMPar.d
            str_par(i,:) = cellstr(strcat('$\;',Par_info.names(i),'$'));
        end
    otherwise % User did not specify names: use x_1, x_2, etc. instead
        for i = 1:AMALGAMPar.d
            str_par(i,:) = cellstr(strcat('$\;x_{',num2str(i),'}$'));
        end
end
% Now create string for Table
str_table = strrep(str_par,'$',''); str_table = strrep(str_table,'\;','');

% Now take first row of min and max of Par_info
Par_info.min = Par_info.min(1,:); Par_info.max = Par_info.max(1,:);
% Rank final parent population and focus on rank 1 solutions
R = AMALGAM_rank(FX,options); idR = (R == 1);
% Now isolate parameter and objective function values of rank 1 solutions
X1 = X(idR,1:AMALGAMPar.d); FX1 = FX(idR,1:AMALGAMPar.m);
% Calculate mean and standard deviation of Pareto solutions
MEAN = mean(X1,1); STD = std(X1,1,1); MED = median(X1,1);

% Determine optimum of each objective function
MAP = nan(AMALGAMPar.m,AMALGAMPar.d);
for ii = 1:AMALGAMPar.m
    % Sort objective functions in ascending order
    [~,ii_sort] = sort(FX(1:AMALGAMPar.N,ii));
    % Find optimum for single criterion
    MAP(ii,1:AMALGAMPar.d) = X(ii_sort(1),1:AMALGAMPar.d);
end

%<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>
%                  TABLE WITH RESULTS OF MEAN AND STD
%<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>

% open an output file with results of formulation
fid = fopen('AMALGAM_output.txt','w');
fprintf(fid,['--------------------------- AMALGAM output file ---------' ...
    '------------------\n']);
fprintf(fid,'==================================================\n');
fprintf(fid,'   x         MEDIAN       MEAN        STD         \n');
fprintf(fid,'--------------------------------------------------\n');
str = str_table; %plotting.str;
% Now print
for j = 1 : AMALGAMPar.d
    if j < 10
        fprintf(fid,'    %s  %+6.3e  %+6.3e  %+6.3e \n', ...
            char(str(j)),MED(j),MEAN(j),STD(j));
    elseif j >= 10 && j < 100
        fprintf(fid,'   %s  %+6.3e  %+6.3e  %+6.3e \n', ...
            char(str(j)),MED(j),MEAN(j),STD(j));
    else
        fprintf(fid,'  %s  %+6.3e  %+6.3e  %+6.3e \n', ...
            char(str(j)),MED(j),MEAN(j),STD(j));
    end
end
fprintf(fid,'==================================================\n');
fprintf(fid,'\n');

%<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>
%                   TABLE WITH RESULTS OF CORRELATION
%<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>

% Calculate AMALGAMPar.d-dimensional parameter correlation matrix
CORR = corrcoef(X1);
evalstr = strcat('fprintf(fid,'' %8.2f');
for k = 1 : AMALGAMPar.d - 1
    evalstr = strcat(evalstr,' %8.2f');
end
evalstr = strcat(evalstr,' %8.2f\n''');
% Now add corr values
evalstr = strcat(evalstr,',');
for k = 1 : AMALGAMPar.d - 1
    evalstr = strcat(evalstr,'corr(',num2str(k),'),');
end
% Now add corr values
evalstr = strcat(evalstr,'corr(',num2str(AMALGAMPar.d),'));');
fprintf(fid,'\n');
% Now print to screen
fprintf(fid,['======================== CORRELATION COEFFICIENTS ======' ...
    '===================\n']);
fprintf(fid,'%9s',[   ]);
for i = 1 : AMALGAMPar.d
    fprintf(fid,'%9s',char(str(i)));
end
fprintf(fid,'\n');
for i = 1 : AMALGAMPar.d
    fprintf(fid,'%9s',char(str(i))); corr = CORR(i,1:AMALGAMPar.d); %#ok
    eval(char(evalstr)); fprintf(fid,'\n');
end
fprintf(fid,['========================================================' ...
    '===================\n']);
fprintf(fid,['------------------------- End AMALGAM output file ------' ...
    '-------------------\n']);

fprintf(fid,'\n');
fclose(fid);

% Now print to screen or not (not on unix/linux)
if ( ispc || ismac ), edit AMALGAM_output.txt, end

str = str_par;

%<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>
%                   PLOT EMPTY FIGURE FOR PDF FILE
%<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>

figure('units','normalized','outerposition',[0 0 1 1], ...
    'name','AMALGAM: Introductory page visual results', ...
    'numbertitle','off');
plot([],[],'ro'); axis([0 1 0 1]); set(gcf,'color','w');
set(gca,'XColor','w','YColor','w');
text(0.3*x_mult,0.6*y_mult,'Visual results of AMALGAM toolbox', ...
    'fontsize',fontsize_titlepage,'interpreter','latex');
% Now about Tables
text(0.3*x_mult,0.5*y_mult,'$\;\;\;\;\;$Tables are not print to PDF file', ...
   'fontsize', fontsize_titlepage,'interpreter','latex');

%<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>
%         PLOT THREE-DIMENSIONAL SNAPSHOT OF PARETO FRONT
%<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>

if ( AMALGAMPar.m == 3 ) % --> only with three objectives
    % open figure
    figure('units','normalized','name',['AMALGAM: Three-dimensional ' ...
        'snapshots of Pareto front'],'numbertitle','off', ...
        'outerposition',[0 0 1 1]);
    % plot pareto solution set
    if ~isempty(Fpareto)
        plot3(Fpareto(:,1),Fpareto(:,2),Fpareto(:,3),'ko', ...
            'markersize',7,'linewidth',2,'MarkerEdgeColor','k', ...
            'MarkerFaceColor','k'); hold on;
        plot3(FX1(:,1),FX1(:,2),FX1(:,3),'ks', ...
            'color',[0.75 0.75 0.75],'MarkerEdgeColor', ...
            [0.75 0.75 0.75],'MarkerFaceColor',[0.75 0.75 0.75]);
    else
        plot3(FX1(:,1),FX1(:,2),FX1(:,3),'ks', ...
            'color',[0.75 0.75 0.75],'MarkerEdgeColor', ...
            [0.75 0.75 0.75],'MarkerFaceColor',[0.75 0.75 0.75]); hold on;
    end
    % scaling of figure
    a = axis; ra = 0.2*(a(2) - a(1)); a_min = a(1) - ra;
    a_max = a(2) + ra; rb = 0.2*(a(4) - a(3)); b_min = a(3) - rb;
    b_max = a(4) + rb;
    rc = 0.2*(a(6) - a(5)); c_min = a(5) - rc; c_max = a(6) + rc;
    axis([a_min a_max b_min b_max c_min c_max]);
    % Find the single criterion ends
    [~,id1] = min(FX1(:,1)); [~,id2] = min(FX1(:,2));
    [~,id3] = min(FX1(:,3));
    % Plot these ends
    plot3(FX1(id1,1),FX1(id1,2),FX1(id1,3),'rx', ...
        'markersize',15,'linewidth',2);
    plot3(FX1(id2,1),FX1(id2,2),FX1(id2,3),'bx', ...
        'markersize',15,'linewidth',2);
    plot3(FX1(id3,1),FX1(id3,2),FX1(id3,3),'gx', ...
        'markersize',15,'linewidth',2);

    % Increase fontsize of numbers
    set(gca,'fontsize',18); grid on
    % Add xlabel
    xlabel('$f_{1}$','fontsize',fontsize_xylabel,'interpreter','latex');
    % Add ylabel
    ylabel('$f_{2}$','fontsize',fontsize_xylabel,'interpreter','latex');
    % Add zlabel
    zlabel('$f_{3}$','fontsize',fontsize_xylabel,'interpreter','latex');
    % Add title
    title(['AMALGAM: Three-dimensional plot of Pareto solution set in ' ...
        'objective function space'],'fontsize',fontsize_title, ...
        'interpreter','latex');
    % And legend strings, optima
    lbl1 = strcat('{\color{red} F_{',num2str(1),',opt}}');
    lbl2 = strcat('{\color{blue} F_{',num2str(2),',opt}}');
    lbl3 = strcat('{\color{green} F_{',num2str(3),',opt}}');
    % Add legend
    if ~isempty(Fpareto)
        legend('{\color{black} $\rm{true}$}',...
            '{\color{gray} AMALGAM}',lbl1,lbl2,lbl3, ...
            'fontsize',fontsize_legend,'box','off');
        %        try set(objh,'interpreter','latex','fontsize',16); catch ME; end
        %        set(legh,'interpreter','latex');
        %set(objh,'linewidth',2);
        %        legh.FontSize = 20;
    else
        legend('{\color{gray} AMALGAM}',lbl1,lbl2,lbl3, ...
            'fontsize',fontsize_legend,'box','off');
        %        try set(objh,'interpreter','latex','fontsize',16); catch ME; end
        %set(objh,'linewidth',2);
        %        set(legh,'interpreter','latex');
        %        legh.FontSize = 20;
    end

end

%<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>
%          PLOT TWO-DIMENSIONAL SNAPSHOTS OF PARETO FRONT
%<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>

% now create all possible combinations
C = nchoosek(1:AMALGAMPar.m,2); n_C = size(C,1); counter = 0;
% Now plot all those
for i = 1:n_C
    % check counter
    if ( counter == 0 ) || ( counter == 6 )
        % open new figure
        figure('units','normalized','name',['AMALGAM: Two-dimensional ' ...
            'snapshots of Pareto front'],'numbertitle','off', ...
            'outerposition',[0 0 1 1]);
        % reset counter
        counter = 0;
    end
    % create subplot
    subplot(2,3,counter + 1);
    % isolate dimensions
    d1 = C(i,1); d2 = C(i,2);
    % plot pareto solution set
    if ~isempty(Fpareto)
        % True nondominated solutions
        plot(Fpareto(:,d1),Fpareto(:,d2),'ko','markersize',7, ...
            'linewidth',2,'MarkerEdgeColor','k', ...
            'MarkerFaceColor','k'); hold on;
        % Rank 1 solutions
        plot(FX1(:,d1),FX1(:,d2),'ks', ...
            'color',[0.75 0.75 0.75], ...
            'MarkerEdgeColor',[0.75 0.75 0.75], ...
            'MarkerFaceColor',[0.75 0.75 0.75]);
    else
        % Rank 1 solutions
        plot(FX1(:,d1),FX1(:,d2),'ks', ...
            'color',[0.75 0.75 0.75], ...
            'MarkerEdgeColor',[0.75 0.75 0.75], ...
            'MarkerFaceColor',[0.75 0.75 0.75]); hold on;
    end
    % scaling of figure
    a = axis; ra = 0.2*(a(2) - a(1)); a_min = a(1) - ra;
    a_max = a(2) + ra; rb = 0.2*(a(4) - a(3));
    b_min = a(3) - rb; b_max = a(4) + rb;
    axis([a_min a_max b_min b_max]);
    % Find the single criterion ends
    [~,id1] = min(FX1(:,d1)); [~,id2] = min(FX1(:,d2));
    % Plot these ends
    evalstr_color_1 = strcat(symbol(d1),'+');
    evalstr_color_2 = strcat(symbol(d2),'+');
    plot(FX1(id1,d1),FX1(id1,d2),char(evalstr_color_1), ...
        'markersize',12,'linewidth',3);
    plot(FX1(id2,d1),FX1(id2,d2),char(evalstr_color_2), ...
        'markersize',12,'linewidth',3);
    % Increase fontsize of numbers
    set(gca,'fontsize',16);
    % Add xlabel
    evalstr = strcat('$f_{',num2str(d1),'}$');
    xlabel(evalstr,'fontsize',fontsize_xylabel,'interpreter','latex');
    % Add ylabel
    evalstr = strcat('$f_{',num2str(d2),'}$');
    ylabel(evalstr,'fontsize',fontsize_xylabel,'interpreter','latex');
    % And legend strings, optima
    lbl1 = strcat('{\color{',char(symbol_text(d1)),['} ' ...
        'F_{'],num2str(d1),',opt}}');
    lbl2 = strcat('{\color{',char(symbol_text(d2)),['} ' ...
        'F_{'],num2str(d2),',opt}}');
    % Add legend if Fpareto exists
    if ~isempty(Fpareto)
        % Add legend
        legend('{\color{black} True front}', ...
            '{\color{gray} Rank 1}',lbl1,lbl2,'box','off', ...
            'location','NorthEast', ...
            'fontsize',fontsize_legend);
        %        set(legh,'interpreter','latex');
        %        set(objh,'linewidth',2);
    else
        % Add legend
        legend('{\color{gray} Rank 1}',lbl1,lbl2, ...
            'box','off','location','NorthEast', ...
            'fontsize',fontsize_legend);
        %        set(legh,'interpreter','latex');
        %        set(objh,'linewidth',2);
    end
    if (counter == 1)
        evalstr = strcat(['AMALGAM: Two-dimensional snapshots of ' ...
            'Pareto samples in objective function space']);
        %shg; %mtit(char(evalstr),'fontsize',18);
        title(char(evalstr),'fontsize',fontsize_title);
    end
    % update counter;
    counter = counter + 1;
end

%<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>
%                PLOT MAP WITH CORRELATION ESTIMATES
%<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>

if ismember(AMALGAMPar.d,2:50)
    figure('units','normalized','name',['AMALGAM: Pareto correlation ' ...
        'matrix'],'numbertitle','off', ...
        'outerposition',[0 0 scr_z(4)/scr_z(3)+0.05 1]);
    % Plot
    ax_corr = axes('units','normalized'); axpos_corr = ...
        [ 0.06 0.06 0.9 0.9 ];
    set(ax_corr,'position',axpos_corr);
    imagesc(ax_corr,CORR); hcol = colorbar; set(hcol,'tickdir','out',...
        'fontsize',fontsize_axis_numbers);
    tickvalues = get(hcol,'Ticks'); tickvalues_fine = cell(1, ...
        numel(tickvalues));
    for i = 1:numel(tickvalues)
        tickvalues_fine{i} = num2str(tickvalues(i),'%3.2f');
    end
    set(hcol,'ticks',tickvalues,'ticklabels',tickvalues_fine,'tickdir', ...
        'out');
    % Now change axis
    set(ax_corr,'fontsize',fontsize_axis_numbers);
    % adjust position of ticks
    set(ax_corr,'XTick',1:AMALGAMPar.d+1,'xticklabel',[]);
    set(ax_corr,'YTick',1:AMALGAMPar.d+1,'yticklabel',[]);
    for j = 1 : AMALGAMPar.d
        % Now add labels as well
        h = text(ax_corr,j,AMALGAMPar.d + .5 + 0.055*AMALGAMPar.d, ...
            str_par(j),'fontsize',fontsize_labels-2);
        set(h, 'rotation', 90)
        text(ax_corr,0.5 -0.055*AMALGAMPar.d,j,str_par(j),'fontsize',...
            fontsize_labels-2);
    end
    % set labels
    set(ax_corr,'xtick', linspace(0.5,AMALGAMPar.d-0.5,AMALGAMPar.d), ...
        'ytick',linspace(0.5,AMALGAMPar.d-0.5,AMALGAMPar.d));
    set(ax_corr,'xgrid', 'on', 'ygrid', 'on', 'gridlinestyle', '-', ...
        'xcolor', 'k', 'ycolor', 'k');
    % Add title
    title(ax_corr,strcat(['AMALGAM: Map of correlation coefficients of ' ...
        'Pareto parameter samples']),'fontsize',fontsize_title);
    % Set font
    set(ax_corr,'fontsize',fontsize_axis_numbers); set(gcf,'color','w');
    % Plot box around figure;
    %plot_box(ax_corr);
elseif ( AMALGAMPar.d == 1 )
    fprintf('\n');
    fprintf(['AMALGAM WARNING: Cannot plot map with bivariate scatter ' ...
        'plots as AMALGAMPAR.d = 1\n']);
else
    fprintf('\n');
    fprintf(['AMALGAM WARNING: Cannot plot map with Pareto corelation ' ...
        'estimates sd AS AMALGAMPAR.d = %1d [= too large]\n'], ...
        AMALGAMPar.d);
end

%<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>
%           CORRELATION PLOTS OF THE PARETO PARAMETER SAMPLES
%<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>

% Only plot this matrix if less or equal to 30 parameters
if ismember(AMALGAMPar.d,1:25)
    % Open a new plot
    figure('units','normalized','name',['AMALGAM: Marginal ' ...
        'distribution and bivariate scatter plots of Pareto samples'], ...
        'numbertitle','off','outerposition',[0 0 1 1]);
    % Plot a matrix (includes unscaled marginals on main diagonal!
    [H,AX,~,P,PAx] = plotmatrix(X1(:,1:AMALGAMPar.d),'+r'); hold on;
    % Now increase font size of each figure - except main diagonal
    set(AX,'fontsize',12);
    % Add title
    title(['AMALGAM: Marginal distribution and bivariate scatter ' ...
        'plots of Pareto samples'],'fontsize',fontsize_title);
    % label the plots
    for i = 1:length(AX)
        ylabel(AX(i,1),str(i),'fontsize',fontsize_xylabel);
        xlabel(AX(end,i),str(i),'fontsize',fontsize_xylabel);
    end
    % Now determine ranges
    %    mm1 = minmax(Pars(:,1:DREAMPar.d)'); dm = mm1(:,2) - mm1(:,1);
    mm1 = [ min(X1(:,1:AMALGAMPar.d))' max(X1(:,1:AMALGAMPar.d))' ];
    dm = mm1(:,2) - mm1(:,1);
    mm(:,1) = mm1(:,1) - dm/10; mm(:,2) = mm1(:,2) + dm/10;
    rr = 10.^(2 - order_no(mm)); mm = round( rr.*mm )./rr;
    minP = mm(:,1); maxP = mm(:,2);
    idx_not = find(minP == maxP); minP(idx_not) = mm1(idx_not,1);
    maxP(idx_not) = mm1(idx_not,2);

    % Now make sure that there are not too many bins in histograms
    for i = 1: AMALGAMPar.d
        if verLessThan('matlab','9.1')
            set(P(i),'Facecolor',[0.5 0.5 0.5],'EdgeColor','w');
        else
            set(P(i),'NumBins',20,'Facecolor',[0.5 0.5 0.5], ...
                'EdgeColor','w','BinLimits',[minP(i) maxP(i)]);
            %,'BinWidth',[maxP(i)-minP(i)]/20);
        end
        set(PAx(i),'Xlim',[minP(i) maxP(i)]);
    end
    for i = 1: AMALGAMPar.d * AMALGAMPar.d
        set(H(i),'Marker','.','markersize',12,'color',[0.5 0.5 0.5]);
    end
    % Now change all axis to lie between minP and maxP
    for i = 1 : AMALGAMPar.d
        for j = 1 : AMALGAMPar.d
            hold(AX(j,i),'on');
            if ( j ~= i )
                if verLessThan('matlab','9.1')
                    % do nothing
                else
                    h = lsline(AX(j,i)); set(h,'color','k', ...
                        'linestyle','--','linewidth',1);
                end
            end
            set(AX(j,i),'Xlim',[minP(i) maxP(i)],'Ylim', ...
                [minP(j) maxP(j)]); %axis tight
        end
    end
    % Now add the prior ranges - if hypercube
    if isfield(Par_info,'boundhandling')
        for i = 1: AMALGAMPar.d
            for j = 1 : AMALGAMPar.d
                if i ~= j
                    hold(AX(j,i),'on');
                    % Vertical lines
                    plot(AX(j,i),[ Par_info.min(i) ...
                        Par_info.min(i) ],[ Par_info.min(j) ...
                        Par_info.max(j)],'b','color',[0 0.4470 0.7410]);
                    plot(AX(j,i),[ Par_info.max(i) ...
                        Par_info.max(i) ],[ Par_info.min(j) ...
                        Par_info.max(j)],'b','color',[0 0.4470 0.7410]);
                    % Horizontal lines
                    plot(AX(j,i),[ Par_info.min(i) ...
                        Par_info.max(i) ],[ Par_info.min(j) ...
                        Par_info.min(j)],'b','color',[0 0.4470 0.7410]);
                    plot(AX(j,i),[ Par_info.min(i) ...
                        Par_info.max(i) ],[ Par_info.max(j) ...
                        Par_info.max(j)],'b','color',[0 0.4470 0.7410]);
                elseif (i == j)
                    % hold(AX(i,j),'on');
                    % AX(i,j); line([ Par_info.min(i), Par_info.min(i)], ...
                    %    [0 1.5*max(P(i).Values)],'linewidth',1, ...
                    %    'color',[0 0.4470 0.7410]);
                    % set(h,'linestyle',':','linewidth',2, ...
                    %    'color',[0.5 0.5 0.5]);
                    % line([ Par_info.max(i), Par_info.max(i)], ...
                    %    ylim,'linewidth',1,'color',[0 0.4470 0.7410]);
                    % set(h,'linestyle',':','linewidth',2, ...
                    %    'color',[0.5 0.5 0.5]);
                end
            end
        end
    end
else
    fprintf('\n');
    fprintf(['AMALGAM WARNING: Cannot plot bivariate scatter plots ' ...
        'as AMALGAMPAR.d = %1d (= too large)\n'],AMALGAMPar.d);
end

%<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>
%                CONVERGENCE TO PARETO DISTRIBUTION
%<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>

if ~isempty(Fpareto)
    % Open new figures
    figure('units','normalized','name',['AMALGAM: Convergence to ' ...
        'Pareto distribution'],'numbertitle', ...
        'off','outerposition',[0 0 1 1]);
    % Plot the evolution of the IGD statistic
    semilogy(output.IGD(:,1),output.IGD(:,2),'r');
    % Increase fontsize
    set(gca,'fontsize',fontsize_axis_numbers);
    % Now play with xticks
    axis([0 output.IGD(end,1) 0.5 * min(output.IGD(:,2)) ...
        1.2*max(output.IGD(:,2))]);
    % Add a title
    title('AMALGAM: Evolution of IGD statistic/diagnostic', ...
        'fontsize',fontsize_title,'interpreter','latex');
    % Add xlabel
    xlabel('Number of generations','fontsize',fontsize_xylabel, ...
        'interpreter','latex');
    % Add Ylabel
    ylabel('Inverted Generational Distance, IGD', ...
        'fontsize',fontsize_xylabel,'interpreter','latex');
end

%<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>
%                 PLOT THE PARETO PARAMETER SOLUTIONS
%<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>

% How many rank 1 solutions do we have?
Nrank = size(X1,1);
% First normalize each parameter with respect to its prior ranges
X1_n = ( X1 - repmat( Par_info.min,Nrank,1) )./ ...
    ( repmat ( Par_info.max - Par_info.min , Nrank,1) );
% Locate the best individual solutions for each objective
r = nan(AMALGAMPar.m,1);
for j = 1 : AMALGAMPar.m
    % Find best solution
    i = find(FX1(:,j) == min(FX1(:,j))); r(j,1) = i(1);
end
figure('units','normalized','name',['AMALGAM: Normalized Pareto ' ...
    'parameter ranges'],'numbertitle','off', ...
    'outerposition',[ 0.05 0.3 0.9 0.6 ]);
% % Open new figure
% figure('units','normalized','outerposition',[0 0 1 1]);

% Now plot each Pareto solution in normalized space
h = plot(1:AMALGAMPar.d,X1_n(:,1:AMALGAMPar.d), ...
    'color',[0.75 0.75 0.75]); hold on;
% Define color
set(h(1),'color',[0.75 0.75 0.75]);
% Define colors for different single criterion optima
symbol = {'r','b','g','c','k','m','y'};
% Now add single criterion "best" solutions;
for j = 1 : AMALGAMPar.m
    % Plot best solution for jth objective
    plot(1:AMALGAMPar.d,X1_n(r(j),1:AMALGAMPar.d), ...
        char(symbol(j)),'linewidth',2);
    % Define color
    set(h(j+1),'color',char(symbol(j)));
end
% create legend
leg_string = cell(1,AMALGAMPar.m+1);
leg_string(1) = {'Rank 1'};
% Add legend
for i = 1:AMALGAMPar.m
    leg_string(i+1) = {[ '$F_{' , num2str(i) , ',\rm{opt} }$']};
end
% Legend
[legh,objh] = legend(leg_string,'box','off','fontsize',20);         %#ok
set(objh,'linewidth',3);
try set(objh,'interpreter','latex','fontsize',16); catch, end
% Now plot vertical lines
for j = 1 : AMALGAMPar.d
    % Plot best solution for jth objective
    plot([j j],[0 1],'k--','linewidth',1);
    % Now add labels as well
    h = text(j,-0.09,str(j),'fontsize',20,'interpreter','latex');
    set(h, 'rotation', 90)
end
% Define axis
axis([0.8 AMALGAMPar.d + 0.2 0 1]);
% Now add parameter labels to X-axis.
set(gca,'xticklabel',[]); set(gca,'fontsize',fontsize_axis_numbers);
% Then add title
title('AMALGAM: Normalized ranges of Pareto solution samples', ...
    'fontsize',fontsize_title,'interpreter','latex');
% Define ylabel
ylabel('Normalized ranges','fontsize',fontsize_xylabel, ...
    'interpreter','latex');

%<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>
%              PLOT THE CONTRIBUTION OF EACH ALGORITHM
%<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>

% Open new figure
figure('units','normalized','name',['AMALGAM: Selection probabilities ' ...
    'of recombination methods'],'numbertitle','off','outerposition', ...
    [ 0.05 0.3 0.9 0.6 ]);
% Plot contribution of each algorithm
ha = plot(output.p_rm(:,1),output.p_rm(:,2:end),'linewidth',2);
% Increase fontsize
set(gca,'fontsize',fontsize_axis_numbers);
% Now change colors
for i = 1:numel(ha)
    set(ha(i),'color',colororder(i,1:3));
end
% Return maximum selection probability to scale y-axis
M = max(max(output.p_rm(:,2:end)));
% Now play with xticks
axis([0 output.p_rm(end,1) 0 min(1,1.1*M)]);
% Add a title
title('AMALGAM: Selection probability of each recombination method', ...
    'fontsize',fontsize_title,'interpreter','latex');
% Add xlabel
xlabel('Number of generations','fontsize',fontsize_xylabel, ...
    'interpreter','latex');
% Add ylabel
ylabel('Selection probability','fontsize',fontsize_xylabel, ...
    'interpreter','latex');
% Legend
legend(upper(AMALGAMPar.rec_methods),'fontsize',fontsize_legend, ...
    'interpreter','latex','location','NorthEast','box','off');
%set(objh,'linewidth',2);

%<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>
%             PLOT THE PARETO SIMULATION UNCERTAINTY
%<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>

if ~isempty(YX)
    % Determine min and max at each time
    sim_min = min(YX); sim_max = max(YX);
    % Open new figure
    figure('units','normalized','name',['AMALGAM: Pareto Simulation ' ...
        'uncertainty'],'numbertitle','off', ...
        'outerposition',[ 0.05 0.3 0.9 0.6 ]);
    % We start with the total uncertainty
    Fill_Ranges(1:size(YX,2),sim_min,sim_max,[0.75 0.75 0.75]); hold on;
    % Determine y-range
    delta_y = max(sim_max) - min(sim_min);
    if min(sim_min) > 10
        % Fit axes
        axis([0 size(YX,2) min(sim_min) - 0.1*delta_y ...
            max(sim_max) + 0.1*delta_y]); set(gca,'fontsize',16);
    else
        % Fit axes
        axis([0 size(YX,2) 0 max(sim_max) + 0.1*delta_y ]);
        set(gca,'fontsize',16);
    end
    % Add xlabel
    xlabel('Time, $t$','fontsize',18,'interpreter','latex');
    % Add ylabel
    ylabel('Variable of interest, $y$','fontsize',18,'interpreter','latex');
    % Add title
    title('AMALGAM: Pareto simulation uncertainty ranges', ...
        'fontsize',fontsize_title,'interpreter','latex');
    % Legend
    [legh] = legend('Pareto simulation uncertainty','fontsize',16, ...
        'location','NorthEast','box','off');
    set(legh,'interpreter','latex');
end

%<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>
%             PLOT THE DISTRIBUTION OF EACH PARAMETER
%<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>

% Determine id of rank 1 solutions of Z (all populations)
id = rank_Z(Z,AMALGAMPar,options);
% Extract parameter values of rank 1 solutions Z
XZ1 = Z(id,1:AMALGAMPar.d);
% Extract objective function values of rank 1 solutions Z
% FZ1 = Z(id,AMALGAMPar.d+1:AMALGAMPar.d+AMALGAMPar.m);

n_row = 2; n_col = 4;
% Now determine row and column number of each parameter
[row_par,col_par,idx_fig] = deal(nan(AMALGAMPar.d,1));
[row_par(1),col_par(1),idx_fig(1)] = deal(1);
for d = 2:AMALGAMPar.d
    col_par(d) = col_par(d-1)+1; row_par(d) = row_par(d-1);
    if col_par(d) == n_col + 1
        row_par(d) = row_par(d-1)+1; col_par(d) = 1;
    end
    if row_par(d) == n_row + 1
        row_par(d) = 1; idx_fig(d) = 1;
    else
        idx_fig(d) = 0;
    end
end

% Compute number of bins based on different rules
Nbins = nan(1,AMALGAMPar.d);
for i = 1:AMALGAMPar.d
    Nbins(i) = calcnbins(XZ1(:,i));
end

% Take the minimum of the number of bins for each parameter
nbins = min(min(Nbins),maxbins);
% Added 2024 as there are too many bins
nbins = max(5,round(nbins/2));
% End added

% Create title
title_str = strcat('AMALGAM: Marginal Pareto distribution of parameters');
% Create legend
leg_lbl = cell(1,AMALGAMPar.m+1); leg_lbl(1) = {'bar'};
for ii = 1:AMALGAMPar.m
    % And legend strings, optima
    leg_lbl(ii+1) = strcat('{\color{',symbol_text(ii),'} F_{', ...
        num2str(ii),',opt}}');
end
% Now loop over parameters
clear ax_hist
for par = 1:AMALGAMPar.d 
    if idx_fig(par) == 1            % Open new figure/add title to previous one
        if exist('ax_hist','var')   % Add title to previous figure
            warning off
            mtit(fig,char(title_str),'fontsize',fontsize_title,...
                'interpreter','latex'); warning on
        end
        fig = figure('units','normalized','name',['AMALGAM: Marginal ' ...
            'distributions of Pareto parameters [= rank 1]'], ...
            'numbertitle','off','outerposition',[0 0 1 1]);
    end
    % Define first axis
    ax_hist = axes('units','normalized'); axpos_hist = [ ...
        0.06 + (col_par(par)-1) * 0.24 , ...
        0.56 - (row_par(par)-1) * 0.48 , 0.20 , 0.38 ];
    set(ax_hist,'position',axpos_hist);
    % New implementation using histcounts
    [M,edges] = histcounts(XZ1(:,par),nbins,'normalization',...
        'countdensity');
    X = 1/2*(edges(1:nbins)+edges(2:nbins+1)); % midpoint of each bin
    % And plot histogram in red
    h = bar(ax_hist,X,M./max(M),'r'); hold on; 
    % --> can be scaled to 1 if using "trapz(X,N)" instead of "sum(N)"!
    set(h,'Facecolor',[0.5 0.5 0.5],'EdgeColor','w');
    % Now determine the ranges of X
    set(ax_hist,'fontsize',fontsize_axis_numbers,'tickdir','out',...
        'ticklength',[0.03 0.05],'XMinorTick','on',...
        'YMinorTick','on'); axis tight;
    % Define xtick values
    xlabel(ax_hist,str_par(par),'fontsize',fontsize_labels,...
        'interpreter','latex');
    % Add legend first subplot
    switch col_par(par)
        case 1  % Add y-label
            set(ax_hist,'ytick',0:0.2:1.0,'yticklabel',...
                {'0.0','0.2','0.4','0.6','0.8','1.0'});
            ylabel(ax_hist,'Empirical density','fontsize',fontsize_labels);
        otherwise
            set(ax_hist,'ytick',0:0.2:1.0,'yticklabel',[]);
    end
    % Lets add the MAP value - each objective
    for ii = 1:AMALGAMPar.m
        plot(ax_hist,MAP(ii,par),0,char(strcat(symbol(ii),'x')), ...
            'Markersize',25,'linewidth',linewidth_marker);
    end
    fix_ticklabels(ax_hist,'x');
    if idx_fig(par) == 1  % print legend
        [~, hico] = legend(ax_hist,leg_lbl,'fontsize',fontsize_legend, ...
            'box','off');
        % Delete objects associated with last 2 black lines
        istxt = strcmp(get(hico, 'type'), 'text');
        hicot = hico(istxt);
        hicol = hico(~istxt);
        delete(hicot(ismember(get(hicot, 'String'), ...
            {'bar','data1','data2'})));
        delete(hicol(ismember(get(hicol, 'Tag'),    ...
            {'bar','data1','data2'})));
    end
    % Adjust the axis
    %    axis(ax_hist,[x_min x_max 0 maxY]);
    % Add label
    fig_code = strcat('(',ranktoletter(par),')');
    text(ax_hist,0.02,0.94,fig_code,'units','normalized',...
        'fontsize',fontsize_A,'interpreter','latex',...
        'horizontalalignment','left');
    % If DREAMPar.d then add title
    if par == AMALGAMPar.d
        mtit(fig,char(title_str),'fontsize',fontsize_title,...
            'interpreter','latex');
    end
    % Plot box around figure;
    plot_box(ax_hist); set(gcf,'color','w');
end

%<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>
%<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>

% Get figure handles
figHandles = flipud(findall(0,'Type','figure'));
for zz = 1:numel(figHandles)
    figure(figHandles(zz)); set(gcf,'color','w');
    switch zz
        case 1
            exportgraphics(figHandles(1),file_name);
        otherwise
            % Append to existing figure
            exportgraphics(figHandles(zz),file_name,'Append',true);
    end
end
% Open PDF document
open(file_name);

end
