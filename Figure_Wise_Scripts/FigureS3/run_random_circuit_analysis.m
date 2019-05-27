%% Set Paths and Load Forward Solutions
% General Parameters
clear; clc; close all;
prefix = 'SM04';                                                            %SM04, cat_004, manny 
meastype = 'meg';                                                           %meg or meg_eeg
add_allpaths;
meas_fname = strcat(prefix,'_SomSens_SingleNoTask_offl_SSP.fif');
chtype = 'all';                                                             %1 for magnetom, 2 for gradiom
gradmag = '1g25m';                                                          %relative ratio gradiom to magnetom
% Forward Files
csvds_fname = strcat(datapath, 'fwd/', prefix, '_cort_meg_', chtype, '_', gradmag,'_SVDs.mat');
ssvds_fname_old = strcat(datapath, 'fwd/', '2014-15_Vectorview_System/whiten_',...
    gradmag, '/',prefix,'_subc_',meastype,'_',chtype,'_', gradmag,'_SVDs.mat'); %2014 Thesis Version
ssvds_fname = strcat(datapath,'fwd/',prefix,'_subc_ico3_',meastype, '_',chtype,'_',gradmag,'_SVDs.mat');
% Parameters of Random Draws Regions to Test
ico_num = 2; num_cort = 5; num_draws = 50;                                  %Cortical Draws Settings
shuffle = 1;                                                                %shuffle random number seeds or not
sregions = {'lvpl','lstr','lhip','ramy','bsred'};                           %scenario name
regname = {'lt','lp','lhipsurf','ra','bsred'};                              %anatomic regions in source space file
regdiv = [1 20 35 7 5];                                                     %divisions of above regions in source space file
type_flag = 'all';                                                          %type of combinations to consider for angles
saveon = 1;                                                                 %to save figures

for scenar = 5:5%1:length(sregions)
if shuffle 
    rng('shuffle'); disp('shuffle');                                        %random number settings
else
    rng('default'); disp('reset');                                          %random number settings                                                        
end                        
slabel_fnames = sregions{scenar};                                           %scenario name
sregname = regname{scenar}; sregsdiv = regdiv(scenar);                      %specifics for loading division of source space 

%% Read in SVD Files
load(csvds_fname,'source','ico','p');                                       %eigendecompositions of cortical patch fwds
cortp = p(ico_num); clear p;                                                %# cortical eigmodes
cpatch = ico(ico_num); num_patch = length(cpatch.patch); clear cpatch;      %cortical patches
if strcmp(slabel_fnames,'lvpl')
    load(ssvds_fname_old,'sdiv_fwd');                                       %eigendecompositions of subcortical volume fwds
else
    load(ssvds_fname,'sdiv_fwd');                                           %eigendecompositions of subcortical volume fwds
end

%% Initialize/Circuit Setup
circ_name = struct('overall',[], 'scort',[],'cort',[]);                     %Initialize Structure Array
circ_name.overall = strcat('random-',slabel_fnames);                        %Random Circuit
circ_name.scort = slabel_fnames;                                            %Name of Subcortical Region
circ_name.cort = cell(1,num_draws);                                         %Random Cortical Patches

%% Subcortical Division
% Pick division location and set up label information
scortgrp = struct('sdivno',[],'sdiv_detail',[], 'sdiv', [], 'left_right', []); %Initialize Structure Arrays
for i = 1:length(sdiv_fwd) 
    all_regions{i} = sdiv_fwd(i).reg_name;                                  % All Region Names
end
ind = find(strcmp(sregname, all_regions));                                  % Identify all region names
if sregsdiv == 0
    sregsdiv = randperm(length(ind),1);                                     % Pick a division at random  
end
slabel_nums = ind(sregsdiv); clear ind;                                     % Identify subdivision index amongst sdiv_fwd divisions
scortgrp.sdivno = slabel_nums;
scortgrp.left_right = strcmp(slabel_fnames(1),'l');                         %1 for left 0 for right
scortgrp.sdiv_detail = sdiv_fwd(scortgrp.sdivno);                           % subcortical reduced decomposition                                      
% Reduced Order Eigendecompositions for Angle Analysis
for ss = 1:length(scortgrp.sdiv_detail)
	sdiv(ss).name = scortgrp.sdiv_detail(ss).reg_name;                      %name of region
    sdiv(ss).hemi = scortgrp.left_right(ss);                                %left hemisphere
    sdiv(ss).sp = scortgrp.sdiv_detail(ss).numwhite_modes;                  %number of modes for NMRA
    sdiv(ss).U = scortgrp.sdiv_detail(ss).whiteUs(:, 1:sdiv(ss).sp);        %left eigenvector
    sdiv(ss).S = scortgrp.sdiv_detail(ss).whiteS(1:sdiv(ss).sp, 1:sdiv(ss).sp); %eigenvalues
    sdiv(ss).V = scortgrp.sdiv_detail(ss).whiteVs(:, 1:sdiv(ss).sp);        %right eigenvector
end
scortgrp.sdiv = sdiv;                                                       %reduced version into group structure

%% Cortical Patches: Random Draws
angs_allcombos_sc = [];     angs_allcombos_cc = [];
cortgrp = struct('patchno', [], 'cpatch_detail', [], 'cpatch', [], 'left_right',[]); %Initialize Structure Arrays

for rdnum = 1:num_draws
disp(strcat('random circuit # ',num2str(rdnum))); tic;  

% Pick patch locations and set up label information
patchno = randperm(num_patch, num_cort);                                    %pick random selection of 4 cortical patches
left_right = logical(patchno <= num_patch/2);                               %1 for left 0 for right
cpatch_detail = ico(ico_num).patch(patchno);                                %cortical reduced decomposition

% Reduced Order Eigendecompositions for Angle Analysis
for cc = 1:length(cpatch_detail)
    cpatch(cc).name = patchno(cc);                                          %patch index
    cpatch(cc).hemi = left_right(cc);                                       %left hemisphere
    cpatch(cc).cp = cortp;                                                  %number of modes for NMRA
    cpatch(cc).U = cpatch_detail(cc).U(:, 1:cpatch(cc).cp);                 %left eigenvector
    cpatch(cc).S = cpatch_detail(cc).S(1:cpatch(cc).cp, 1:cpatch(cc).cp);   %eigenvalues
    cpatch(cc).V = cpatch_detail(cc).V(:, 1:cpatch(cc).cp);                 %right eigenvector
end

% Group Structure
cortgrp(rdnum).patchno = patchno;                                           %store patch numbers considered in circ_name structure
circ_name.cort{rdnum} = patchno;                                            %store patch numbers considered in circ_name structure
cortgrp(rdnum).left_right = left_right;                                     %1 for left 0 for right
cortgrp(rdnum).cpatch_detail = cpatch_detail;                               % cortical reduced decomposition
cortgrp(rdnum).cpatch = cpatch;                                             %reduced version into group structure

%% Compute Angles for All Combinations
%[angs_sc, angs_cc] = compute_angs_allcombos(scortgrp, cortgrp(rdnum), type_flag);
disp('subcort vs cort'); 
angs_sc = reg1_reg2_allcombos(sdiv, cpatch, type_flag, sdiv.sp, cortp);     %angles for all combinations of modes
disp('cort vs cort'); 
angs_cc = reg1_reg2_allcombos(cpatch, cpatch, type_flag, cortp, cortp);     %angles for all combinations of modes
angs_allcombos_sc = [angs_allcombos_sc angs_sc];
angs_allcombos_cc = [angs_allcombos_cc angs_cc];
toc;
clear cpatch_detail cpatch angs_cc angs_sc; 
end
clear sdiv;

%% Plot Angles for All Combinations and Save Results
disp('Source Space Plot');
disp(strcat(sdiv_fwd(scortgrp.sdivno).reg_name, '_submask',num2str(sdiv_fwd(scortgrp.sdivno).index)));
savepath = strcat(pwd, '/', circ_name.overall, '_circ_results/'); 
if shuffle 
savename = strcat(savepath, prefix, '_', circ_name.overall,'_', num2str(num_cort),'c_angles_shuffle');
figsavename = strcat(savepath, circ_name.overall,'_',num2str(num_cort), 'c_allmodecombos_shuffle');
else
savename = strcat(savepath, prefix, '_', circ_name.overall,'_', num2str(num_cort),'c_angles');
figsavename = strcat(savepath, circ_name.overall,'_',num2str(num_cort), 'c_allmodecombos');
end
[summ_stat_cc, summ_stat_sc] = figure_rdraws_results(angs_allcombos_sc, angs_allcombos_cc, saveon, figsavename);
if saveon
    save(savename,'cortgrp','scortgrp','circ_name','angs_allcombos_sc','angs_allcombos_cc', 'chtype', ...
    'ssvds_fname', 'sregname','sregsdiv',...
    'csvds_fname', 'ico_num','num_cort','num_draws');
end

clear slabel_nums sdiv patchno cpatch angs_allcombos_sc angs_allcombos_cc circ_name scortgrp cortgrp sregname sregsdiv
end