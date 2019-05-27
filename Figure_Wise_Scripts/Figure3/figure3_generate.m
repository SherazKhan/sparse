%% Set Paths and Parameters
%cortlabel_extract;
%run_somerp_circuit_analysis;
clear; clc; close all;
addpath(genpath('/autofs/eris/purdongp/users/pavitrak/sourceloc_scripts/ScSPIGH/'));
prefix = 'SM04'; datatype = 'Illustrations'; meastype = 'meg'; add_allpaths;
circ_choice = 'som_erp'; disp(strcat('circuit # ',circ_choice));            %'som_erp', 'random' or 'inf_hip'; %specify circuit 
savepath = strcat(pwd, '/', circ_choice, '_circ_results/'); 
grp_fname = strcat(savepath, 'SM04_forgrpangles_',circ_choice,'_meg_all_1g25m');
load(strcat(grp_fname,'.mat'), 'scortgrp','cortgrp','circ_name');
scortgrp = scortgrp(logical(strcmp({circ_name.overall}, circ_choice)));     %choose the regions corresponding to circuit chosen
cortgrp = cortgrp(logical(strcmp({circ_name.overall}, circ_choice)));       %choose the regions corresponding to circuit chosen
saveon = 1;                                                                 %to save figures
type_flag = 'all';                                                          %type of combinations to consider

%% Read Subcortical and Cortical Components of Circuit and Regions to Analyze
sdiv = scortgrp.sdiv;                                                       %subcortical divisions in circuit
s = length(sdiv); sreg_inds = 1:1:s;                                        %no. and indices of subcortical divisions
[sdiv_subs, sflag] = all_combs(sdiv);                                       %all subsets of subcortical divisions
scort_p = sdiv.sp;
cpatch = cortgrp.cpatch;                                                    %cortical patches in circuit        
[cdiv_subs, cflag] = all_combs(cpatch);                                     %all subsets of subcortical divisions
cort_p = cpatch.cp;

%% Angles for All Combinations of Modes in Source Space: Scort vs. Cort
disp('subcort vs cort'); 
angs_allcombos_sc = reg1_reg2_allcombos(sdiv, cpatch, 'all',scort_p, cort_p);%angles for all combinations of modes
summ_stat_sc = median(angs_allcombos_sc);                                   %summary statistic: median

%% Angles for All Combinations of Modes in Source Space: Cort vs. Cort Control
disp('cort vs cort'); 
angs_allcombos_cc = reg1_reg2_allcombos(cpatch, cpatch, 'all',cort_p, cort_p); %angles for all combinations of modes
summ_stat_cc = median(angs_allcombos_cc);                                   %summary statistic: median

%% Plot and Save Subcort vs. Cort and Cort vs. Cort Histograms
figure, set(gcf,'color','white'); hold all;
%Plot Histograms
num_bins = 20; bw = 90/num_bins;
[~, binlocs_c] = hist(angs_allcombos_cc,num_bins); 
d = diff(binlocs_c)/2;  edges_c = [binlocs_c(1)-d(1), binlocs_c(1:end-1)+d, binlocs_c(end)+d(end)]; edges_c(2:end) = edges_c(2:end)+eps(edges_c(2:end));
[~, binlocs_s] = hist(angs_allcombos_sc,num_bins); 
d = diff(binlocs_s)/2;  edges_s = [binlocs_s(1)-d(1), binlocs_s(1:end-1)+d, binlocs_s(end)+d(end)]; edges_s(2:end) = edges_s(2:end)+eps(edges_s(2:end));
histogram(angs_allcombos_cc, edges_c, 'Normalization','probability', 'EdgeColor','g','FaceColor','g','FaceAlpha', 0.30,'DisplayStyle','bar');%num_bins, 'BinWidth', bw,
histogram(angs_allcombos_sc, edges_s, 'Normalization','probability', 'EdgeColor',[1 0.5 0],'FaceColor',[1 0.5 0],'FaceAlpha', 0.65,'DisplayStyle','bar');%num_bins, 'BinWidth', bw, 
xlim([0 90]); 
% Plot Median Lines
nn = axis; min_no = nn(3); max_no = nn(4);
line(summ_stat_sc*ones(1,30), linspace(min_no, max_no, 30),'LineWidth',3,'color',[1 0.5 0]);
line(summ_stat_cc*ones(1,30), linspace(min_no, max_no, 30),'LineWidth',3,'color','g');
% Save
xlabel('Principal Angle (Degrees)','FontSize',14); ylabel('Normalized Histogram','FontSize',14); 
titlestr{1} = ['SC: Median All Combos All Modes:  ', num2str(summ_stat_sc)]; %display summary in histogram
titlestr{2} = ['CC: Median All Combos All Modes:  ', num2str(summ_stat_cc)]; %display summary in histogram
title(titlestr,'FontSize',14); set(gca,'FontSize',14); xlim([0 90]); 
savename = strcat(savepath, circ_choice,'sc_cc_allcombos_allmodes');
if saveon
    print(gcf,'-depsc',savename);
end
pause;

%%Plot and Save Field Maps for Illustration
fieldmap_saves_allmode;