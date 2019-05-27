function [summ_stat_cc, summ_stat_sc] = figure_rdraws_results(angs_allcombos_sc, angs_allcombos_cc, saveon, savename)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Written by Pavitra Krishnaswamy, May 16 2016
% Save Subcort vs. Cort and Cort vs. Cort Histograms
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Setup Figure Parameters
figure, set(gcf,'color','white'); hold all;
num_bins = 20; bw = 90/num_bins;

%%Plot C vs. C Histogram
[~, binlocs_c] = hist(angs_allcombos_cc,num_bins); 
d = diff(binlocs_c)/2;  edges_c = [binlocs_c(1)-d(1), binlocs_c(1:end-1)+d, binlocs_c(end)+d(end)]; edges_c(2:end) = edges_c(2:end)+eps(edges_c(2:end));
histogram(angs_allcombos_cc, edges_c, 'Normalization','probability', 'EdgeColor','g','FaceColor','g','FaceAlpha', 0.30,'DisplayStyle','bar');%num_bins, 'BinWidth', bw,
xlim([0 90]); 

%%Plot S vs. C Histogram
[~, binlocs_s] = hist(angs_allcombos_sc,num_bins); 
d = diff(binlocs_s)/2;  edges_s = [binlocs_s(1)-d(1), binlocs_s(1:end-1)+d, binlocs_s(end)+d(end)]; edges_s(2:end) = edges_s(2:end)+eps(edges_s(2:end));
histogram(angs_allcombos_sc, edges_s, 'Normalization','probability', 'EdgeColor',[1 0.5 0],'FaceColor',[1 0.5 0],'FaceAlpha', 0.65,'DisplayStyle','bar');%num_bins, 'BinWidth', bw, 
xlim([0 90]); 

%% Compute and Plot Summary Statistics
summ_stat_sc = median(angs_allcombos_sc);                                   %summary statistic: median
summ_stat_cc = median(angs_allcombos_cc);                                   %summary statistic: median
nn = axis; min_no = nn(3); max_no = nn(4);
line(summ_stat_sc*ones(1,30), linspace(min_no, max_no, 30),'LineWidth',3,'color','y');
line(summ_stat_cc*ones(1,30), linspace(min_no, max_no, 30),'LineWidth',3,'color','g');
titlestr{1} = ['SC: Median All Combos All Modes:  ', num2str(summ_stat_sc)]; %display summary in histogram
titlestr{2} = ['CC: Median All Combos All Modes:  ', num2str(summ_stat_cc)]; %display summary in histogram

%% Save Figure
xlabel('Principal Angle (Degrees)','FontSize',14); ylabel('Normalized Histogram','FontSize',14); 
title(titlestr,'FontSize',14); set(gca,'FontSize',14); xlim([0 90]); 
if saveon
    print (gcf,'-dpdf',savename);
    %print(gcf,'-depsc',savename);
end

end