%% Set Paths and Load Forward Solutions
clear; clc; close all;
prefix = 'SM04';                                                            %SM04, cat_004, manny 
meastype = 'meg';                                                           %meg or meg_eeg
add_allpaths;
meas_fname = strcat(prefix,'_SomSens_SingleNoTask_offl_SSP.fif');
chtype = 'all';                                                             %1 for magnetom, 2 for gradiom
gradmag = '1g25m';                                                          %relative ratio gradiom to magnetom
csvds_fname = strcat(datapath, 'fwd/', prefix, '_cort_meg_', chtype, '_', gradmag,'_SVDs.mat');
ssvds_fname = strcat(datapath, 'fwd/', '2014-15_Vectorview_System/whiten_',...
    gradmag, '/',prefix,'_subc_',meastype,'_',chtype,'_', gradmag,'_SVDs.mat'); %2014 Thesis Version
circ_choice = 'som_erp';
savepath = strcat(pwd, '/', circ_choice, '_circ_results/'); 
                                                                   
%% Read in SVD Files
load(csvds_fname,'source','ico','p');                                       %eigendecompositions of cortical patch fwds
ico_num = 2; cortp = p(ico_num); clear p;                                   %# cortical eigmodes
load(ssvds_fname,'sdiv_fwd');                                               %eigendecompositions of subcortical volume fwds

%% Load Relevant Cortical and Subcortical Regions
% Initialize
load(strcat(savepath, 'cregion_grps'),'label_fnames','left_right','lh_patchno','rh_patchno'); %Patch#s Extracted
clabel_fnames = label_fnames; clear label_fnames;                           %Patch#s Extracted
slabel_fnames = {'lvpl'};                                                   %Sdiv Names Extracted
sregname = 'lt'; sregsdiv = 1;                                              % LT Submask 1
for i = 1:length(sdiv_fwd) 
    all_regions{i} = sdiv_fwd(i).reg_name;                                  % All Region Names
end
ind = find(strcmp(sregname, all_regions));                                  % Identify all region names
slabel_nums = ind(sregsdiv); clear ind;                                     % Identify subdivision index amongst sdiv_fwd divisions
circ_name = struct('overall',[],'cort',[],'scort',[]);                      %Initialize Structure Arrays
cortgrp = struct('patchno', [], 'cpatch_detail', [], 'cpatch', [], 'left_right',[]); %Initialize Structure Arrays
scortgrp = struct('sdivno',[],'sdiv_detail',[], 'sdiv', [], 'left_right', []); %Initialize Structure Arrays

% Group Regions into Somatosensory Circuit
circ_name.overall = circ_choice;   
circ_name.cort = clabel_fnames(1:2);     
cortgrp.patchno = [lh_patchno{1}' rh_patchno{2}']';     
cortgrp.left_right = [ones(1, length(lh_patchno{1})) zeros(1, length(rh_patchno{2}))];
circ_name.scort = slabel_fnames;
scortgrp.sdivno = slabel_nums;
scortgrp.left_right = ones(1,length(slabel_nums));            

% Load Relevant Reduced Order Eigendecompositions
cortgrp.cpatch_detail = ico(ico_num).patch(cortgrp.patchno);
scortgrp.sdiv_detail = sdiv_fwd(scortgrp.sdivno);

%% Organize Reduced Order Eigendecompositions for Angle Analysis
% Cortical Eigendecompositions
cpatch_detail = cortgrp.cpatch_detail;
for cc = 1:length(cpatch_detail)
	cpatch(cc).name = cortgrp.patchno(cc);                                  %patch index
    cpatch(cc).hemi = cortgrp.left_right(cc);                               %left hemisphere
    cpatch(cc).cp = cortp;                                                  %number of modes for NMRA
    cpatch(cc).U = cpatch_detail(cc).U(:, 1:cpatch(cc).cp);                 %left eigenvector
    cpatch(cc).S = cpatch_detail(cc).S(1:cpatch(cc).cp, 1:cpatch(cc).cp);   %eigenvalues
    cpatch(cc).V = cpatch_detail(cc).V(:, 1:cpatch(cc).cp);                 %right eigenvector
end
cortgrp.cpatch = cpatch;                                                    %reduced version into group structure

% Subcortical Eigendecompositions
sdiv_detail = scortgrp.sdiv_detail;
for ss = 1:length(sdiv_detail)
	sdiv(ss).name = sdiv_detail(ss).reg_name;                               %name of region
    sdiv(ss).hemi = scortgrp.left_right(ss);                                %left hemisphere
    sdiv(ss).sp = sdiv_detail(ss).numwhite_modes;                           %number of modes for NMRA
    sdiv(ss).U = sdiv_detail(ss).whiteUs(:, 1:sdiv(ss).sp);                 %left eigenvector
    sdiv(ss).S = sdiv_detail(ss).whiteS(1:sdiv(ss).sp, 1:sdiv(ss).sp);      %eigenvalues
    sdiv(ss).V = sdiv_detail(ss).whiteVs(:, 1:sdiv(ss).sp);                 %right eigenvector
end
scortgrp.sdiv = sdiv;
clear cpatch_detail cpatch sdiv_detail sdiv

%% Compute and Plot Angles for Groups
savename = strcat(savepath, prefix, '_forgrpangles_', circ_name.overall,'_', meastype,'_',chtype, '_1g25m.mat');
save(savename,'cortgrp','scortgrp','circ_name', 'chtype', 'csvds_fname', 'ssvds_fname', 'ico_num');
save_plots = 0;                                                             %don't save plots
allcombos_grphist;                                                          %angles: all combinations of modes and divs
%fieldmap_saves;                                                            %sample field maps