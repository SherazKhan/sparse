function reg_matrix_plots(reg_matrix, regname, st, en, saveon, savename, cax, flag_on, cmap)

%% Main Figure
figure, set(gcf,'color','white'); 
imagesc(abs(reg_matrix)); colormap(cmap);
%load('Fig4_colormap'); colormap(Fig4_colormap); brighten(-0.3)
%colormap(flipud(colormap('hot')))
if isempty(cax)
    colorbar;
else
    caxis(cax); colorbar; %axis xy;
end

%% Axis Limits
Xt = st; 	Xl = [1 en(end)]; set(gca, 'XTick', Xt, 'TickDir', 'out', 'Xlim',Xl);
Yt = st; 	Yl = [1 en(end)]; set(gca, 'YTick', Yt, 'TickDir', 'out', 'Ylim',Yl);

%% Settings
if flag_on
    xrotvalue = 45; yrotvalue = 45;                                          % slant
else
    savename = strcat(savename,num2str(flag_on));
    xrotvalue = 0; yrotvalue = 0;                                             % 0,0 if flat OR 45,45 if slant
    grid on;   set(gca,'GridLineStyle','--','GridColor','m','Layer', 'top');    % make it easy to discern regions
    set(gca,'tickdir','out'); set(gca,'ticklength',10*get(gca,'ticklength'));   % Ticks for Placing Text Labels
end

%% Names of Regions
regions = str2mat(regname);
title('Resolution Matrix');
tx = text(Xt, (Yl(end)+10)*ones(1,length(Xt)), regions);                    % Text locations for x-axis
set(tx, 'HorizontalAlignment','center','VerticalAlignment','top','rotation',xrotvalue,'FontSize',14);
set(gca,'XTickLabel','');                                                   % Remove the default labels
ty = text(-10*Xt(1)*ones(1,length(Yt)),Yt+0.5,regions);                     % Text locations for y-axis
set(ty, 'VerticalAlignment','middle','HorizontalAlignment','right','rotation',yrotvalue,'FontSize',14);
set(gca,'YTickLabel','');                                                   % Remove the default labels
clear Xt Yt Xl Yl regions xrotvalue yrotvalue tx ty;        

%% Save Figure with All Correct Rendering Properties
if saveon
    print('-depsc', savename);
    print('-dpdf',savename);
end

end