function hist_figure(choose_ang, saveon, savename, X, ctrl_angle, titlestr)

figure, set(gcf,'color','white');
if numel(X) > 1
    N = hist(choose_ang, X); bar(X, N./sum(N));
else
    [~, X] = hist(choose_ang, X); histnorm(choose_ang, X); 
end
nn = axis; min_no = nn(3); max_no = nn(4);
hold on; line(ctrl_angle*ones(1,30), linspace(min_no, max_no, 30),'LineWidth',3,'color','k');
xlabel('Principal Angle (Degrees)','FontSize',14); ylabel('Normalized Histogram','FontSize',14); 
title(titlestr,'FontSize',14); set(gca,'FontSize',14); xlim([0 90]); 
if saveon
    print(gcf,'-depsc',savename);
end

end