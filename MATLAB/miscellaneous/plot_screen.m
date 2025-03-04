function [fig,ax0] = plot_screen(fig,ax0,AMALGAMPar,F,t,F_par,output)
% Plot to screen the AMALGAM results during each trial

names = {'GA','PSO','AMS','DE'};        % Name of individual methods
symbol = {'rs','gd','mp','co','bh'};    % Define colors for different single criterion optima
%colors = {'r','g','m','c','b'};         % Colors for methods
Colors = [ 1 0 0; 0 0.5 0; 0.5 0.5 0.5; 1 0.64 0; 0 0 1];
gif = 1;                                % Print to screen

if numel(AMALGAMPar.rec_methods) == 1
    switch char(AMALGAMPar.rec_methods)
        case {'ga'}
            idx = 1;
        case {'pso'}
            idx = 2;
        case {'ams'}
            idx = 3;
        case {'de'}
            idx = 4;
    end
    method = char(AMALGAMPar.rec_methods(1));
else
    idx = 5; method = 'AMALGAM';
end

if ( t == 2 )
    [ fig,ax0 ] = plot_ZDT(F_par);
    if gif
        % Capture the plot as an image
        frame = getframe(fig); im = frame2im(frame); [imind,cm] = rgb2ind(im,256);
        imwrite(imind,cm,strcat(method,'_ZDT.gif'),'gif', 'Loopcount',inf);
    end
end
axes(ax0);

h1 = semilogy(F(:,1),F(:,2),char(symbol(idx)),'color',Colors(idx,1:3),'linewidth',1,'markersize',8); axis([0 1 10^-4 10^3]);
h1(1).MarkerFaceColor = Colors(idx,1:3);
IGD_print = num2str(round(10000*output.IGD(t,2))/10000);
N_IGD = numel(IGD_print);
if N_IGD < 6
    IGD_print = strcat(IGD_print,num2str(zeros(1,6-N_IGD)));
end
evalstr_1 = strcat('$ = ',{' '},IGD_print,'$');
h2 = text(0.85,1600,evalstr_1,'interpreter','latex','fontsize',16);

evalstr_2 = strcat('$ = ',{' '},num2str(t),'$');
h0 = text(0.65,1600,evalstr_2,'interpreter','latex','fontsize',16);

axpos1 = [3.0 8.6 1 0.2]; ax1 = axes('units','inches'); set(ax1,'position',axpos1);
h3 = barh(1,t/AMALGAMPar.T); colormap([0.25 0.25 0.25]); axis([0 1 0.5 1.5]); set(gca,'xticklabel',[]); set(gca,'yticklabel',[]);

p_alg = output.p_alg(t,2:AMALGAMPar.q+1); 
if idx < 5 % Individual method was selected
    P_alg = zeros(1,4); P_alg(idx) = p_alg; 
else
    P_alg = p_alg;
end

% If AMALGAM is running with four methods
for qq = 1:4
    axpos1 = [2.6 3.2-(qq-1)*0.35 2.5 0.3]; ax1 = axes('units','inches'); set(ax1,'position',axpos1);
    h4 = barh(qq,P_alg(1,qq)); h4(1).FaceColor = Colors(qq,1:3); %pause
    set(gca,'yticklabel',[]); axis([0 1 -0.5+qq 0.5 + qq]); box on
    if qq < 4
        set(gca,'xticklabel',[]);
    end
    x_loc = [0 0.2 0.40 0.6 0.8 1.0]; x_loc_pr = {'0.0','0.2','0.4','0.6','0.8','1.0'};
    xticks(x_loc);
    if qq == 4
        set(gca,'XTickLabel',[]);
        for k = 1:6
            text(x_loc(k),3.03,x_loc_pr(k),'HorizontalAlignment','Center','fontsize',14) %// Play with the 'ylim(1) -0.1' to place the label as you wish.
        end
        text(-0.023,2.05,'$\mbox{SELECTION PROBABILITY}$','interpreter','latex','fontsize',14);
    end
    h5 = text(-0.26,qq-0.1,upper(names(qq)),'color',Colors(qq,1:3),'rotation',0,'Fontweight','bold','fontsize',14);
end

if gif
    % Capture the plot as an image
    frame = getframe(fig); im = frame2im(frame); [imind,cm] = rgb2ind(im,256);
    imwrite(imind,cm,strcat(method,'_ZDT.gif'),'gif','Writemode','append');
end

delete(h1); delete(h0); delete(h2); delete(h3); delete(h4);% delete(h5);

% Secundairy function
function [ fig,ax0 ] = plot_ZDT(F_par)

% Determine figure size and location + add true Pareto front
fig = figure('units','normalized','name','ZDT function','numbertitle','off','outerposition',[ 0.05 0.05 0.5 0.9 ]);

if ~isempty ( F_par )
    semilogy(F_par(:,1),F_par(:,2),'ko','markersize',2,'linewidth',2,'MarkerEdgeColor','k','MarkerFaceColor','k'); hold on;
end
% scaling of figure
axis square; 
axis([0 1 10^-4 10^3])

pos = get(gca, 'Position');
offset = 0.01;
set(gca, ...
    'Box'         , 'on'                        , ...
    'TickDir'     , 'in'                         , ...
    'XMinorTick'  , 'on'                        , ...
    'YMinorTick'  , 'on'                        , ...
    'TickLength'  , [.02 .02]                    , ...
    'LineWidth'   , 1                            , ...
    'FontSize'    , 17                           , ...
    'Position'    , pos + [0, offset, 0, -offset]);
hold on;
text(0.46,10^-4.6,'$F_{1}(\theta)$','interpreter','latex','fontsize',18); h = text(-0.11,10^-0.8,'$F_{2}(\theta)$','interpreter','latex','fontsize',18);
set(h,'rotation',90); set(gcf, 'color', 'w'); ax0 = gca;

text(0.78,1600,'${\rm IGD}$','interpreter','latex','fontsize',16);
text(0.43,1600,'${\rm GENERATION}$','interpreter','latex','fontsize',16);
text(0.02,1600,'$\mbox{PROGRESS}$','interpreter','latex','fontsize',16);