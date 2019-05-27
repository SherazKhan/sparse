%% Set Paths and Parameters
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

for analtype = 1:3                                                          %for each analysis type considered
  %% Read Subcortical and Cortical Components of Circuit and Regions to Analyze
  sdiv = scortgrp.sdiv;                                                     %subcortical divisions in circuit
  s = length(sdiv); sreg_inds = 1:1:s;                                      %no. and indices of subcortical divisions
  [sdiv_subs, sflag] = all_combs(sdiv);                                     %all subsets of subcortical divisions
  scort_p = sdiv.sp;
  cpatch = cortgrp.cpatch;                                                  %cortical patches in circuit        
  [cdiv_subs, cflag] = all_combs(cpatch);                                   %all subsets of subcortical divisions
  cort_p = cpatch.cp;
  
  %% Set Type of Analysis
  switch analtype
    case 1
        disp('cort vs cort'); savename = 'cc';  no_use_flag = cflag; 
        region1 = cpatch;   region2 = cpatch;   savename1 = '_cc_modes';    %set regions and savename for mode wise analysis
        reg1 = cdiv_subs;   reg2 = cdiv_subs;   savename2 = '_cc_divs';     %set regions and savename for division wise analysis    
        reg1p = cort_p;     reg2p = cort_p; 
    case 2
        disp('subcort vs cort'); savename = 'sc'; no_use_flag = [];
        region1 = sdiv;     region2 = cpatch;   savename1 = '_sc_modes';    %set regions and savename for mode wise analysis
        reg1 = sdiv_subs;   reg2 = cdiv_subs;   savename2 = '_sc_divs';     %set regions and savename for division wise analysis
		reg1p = scort_p;    reg2p = cort_p; 
    case 3
        disp('subcort vs subcort'); savename = 'ss'; no_use_flag = sflag; 
        region1 = sdiv;     region2 = sdiv;     savename1 = '_ss';          %set regions and savename for mode wise analysis
        reg1 = sdiv_subs;   reg2 = sdiv_subs;   savename2 = '_ss_divs';     %set regions and savename for division wise analysis
        reg1p = scort_p;    reg2p = scort_p;
  end
  X = 20; %X = 0:5:90;                                                      %for histogram

  %% Mode Wise Analyses
  % Plot Angles for All Combinations of Modes in Source Space
  angs_allcombos = reg1_reg2_allcombos(region1, region2, 'all',reg1p, reg2p); %angles for all combinations of modes
  summ_stat(1) = median(angs_allcombos);                                    %summary statistic: median
  titlestr{1} = ['Median All Combos All Modes:  ', num2str(summ_stat(1))];  %display summary in histogram
  hist_figure(angs_allcombos, saveon, strcat(savepath, circ_choice, savename1,'_allcombos_allmodes'), X, summ_stat(1), titlestr{1});

  %% Mode Wise Analyses - Only Combinations with Top Mode
  % Plot Angles for All Combinations of Modes in Source Space
  angs_allcombos_top = reg1_reg2_allcombos(region1, region2, 'all_top',reg1p, reg2p); %angles for all combinations of modes
  summ_stat(2) = median(angs_allcombos_top);                                %summary statistic: median
  titlestr{2} = ['Median All Combos All Modes with Top:  ', num2str(summ_stat(2))];  %display summary in histogram
  hist_figure(angs_allcombos_top, saveon, strcat(savepath, circ_choice, savename1,'_allcombos_allmodes_top'), X, summ_stat(2), titlestr{2});

  %% Division Wise Analyses
  % Plot Angles for All Combinations for 1 Randomly Selected Mode Per Region
  angs_allcombos_randmodes = reg1_reg2_allcombos(region1, region2, 'rand',reg1p, reg2p); %angles for all combinations of randomly selected modes
  summ_stat(3) = median(angs_allcombos_randmodes);                          %summary statistic: median
  titlestr{3} = ['Median All Combos Rand Modes:  ', num2str(summ_stat(3))]; %display summary in histogram
  hist_figure(angs_allcombos_randmodes, saveon, strcat(savepath, circ_choice, savename1,'_allcombos_randmodes'), X, summ_stat(3), titlestr{3});

  % Plot Angles for All Combinations for Top-Most Mode Per Region
  angs_allcombos_topmodes = reg1_reg2_allcombos(region1, region2, 'top',reg1p, reg2p);   %angles for all combinations of top modes
  summ_stat(4) = median(angs_allcombos_topmodes);                           %summary statistic: median
  titlestr{4} = ['Median All Combos Top Modes:  ', num2str(summ_stat(4))];  %display summary in histogram
  hist_figure(angs_allcombos_topmodes, saveon, strcat(savepath, circ_choice, savename1,'_allcombos_topmodes'), X, summ_stat(4),titlestr{4});

  % Plot Angles for All Combinations for All Modes per Region (Each Region Represented by All Modes)
  angs_all_combos_divs = [];
  if isequal(analtype, 2)
    for i = 1:length(reg1)
      for j = 1:length(reg2)
        angs_allregs{i,j} = reg1_reg2_all(reg1{i}, reg2{j}, 'all');         %angles for each combination of divisions
        angs_all_combos_divs = [angs_all_combos_divs angs_allregs{i,j}(:)'];%collate angles across all combinations
      end
    end
  else
   for i = 1:length(reg1)
     for j = i+1:length(reg2)
        if no_use_flag(i,j) == 0
        angs_allregs{i,j} = reg1_reg2_all(reg1{i}, reg2{j}, 'all');         %angles for each combination of divisions
        angs_all_combos_divs = [angs_all_combos_divs angs_allregs{i,j}(:)'];%collate angles across all combinations    
        end
     end
   end
  end
  summ_stat(5) = median(angs_all_combos_divs);                              %summary statistic
  titlestr{5} = ['Median All Modes:  ', num2str(summ_stat(5))];             %display summary in histogram
  hist_figure(angs_all_combos_divs, saveon, strcat(savepath, circ_choice, savename2,'_allcombos_allmodes'), X, summ_stat(5), titlestr{5});
  
  %% Save and Clear
  save(strcat(grp_fname, '_', savename, '.mat'),'scortgrp','cortgrp',...
    'angs_allcombos','angs_allcombos_top','angs_allcombos_randmodes','angs_allcombos_topmodes','angs_all_combos_divs','angs_allregs',...
    'X', 'region1', 'region2', 'reg1', 'reg2', 'reg1p', 'reg2p', 'no_use_flag','-v7.3');
  pause; close all; 
end