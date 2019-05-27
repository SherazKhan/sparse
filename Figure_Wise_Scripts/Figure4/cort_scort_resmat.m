%% Set Parameters and Paths
clear; clc; close all;
prefix = 'SM04'; meastype = 'meg'; datatype = 'Illustrations';
g = 1; m = 25; gradmag = strcat(num2str(g), 'g', num2str(m), 'm'); chtype = 'all'; %scaling ratios for gradiom vs. magnetom                                                            
dpath = '/autofs/eris/purdongp/users/pavitrak/'; add_allpaths;              %add all paths
fwd_path = strcat(datapath, '/fwd/');                                       %latest forward path

%% Set Parameters for Analysis
%choose_cort = 'all';        choose_scort = 'all';                          %full cortical and subcortical
choose_cort = 'sparse';    choose_scort = 'all';                            %sparse cortical and subcortical (can alter to rand)
%choose_cort = 'all';       choose_scort = 'none';                          %cortical only
%choose_cort = 'none';      choose_scort = 'all';                           %subcortical only
nummodes = 1;                                                               %#eigenmodes - 1 for one mode, 0 for all modes per division
SNR = 9;                                                                    %standard power SNR for ERPs - vary 9 to 25
N = 306;                                                                    %# channels being used
cdiv = 2;                                                                   %the subdiv to consider for this stage

%% Load Whitened Cortical SVDs
csvds_fname = strcat(fwd_path, prefix, '_cort_meg_', chtype, '_', gradmag, '_SVDs.mat');
load(csvds_fname, 'ico','p');                                               %results of cortical patch decomposition
ico = ico(cdiv);                                                            %the ico we need to consider
cp = p(cdiv); clear p                                                       %# eigenmodes for cortical stage
load(strcat(fwd_path,'curr_den.mat'),'cort_cd');                            %[Am/mm^2] - cortical current density
nc = length(ico.patch);                                                     %for ico-cdiv
lk = select_cortpatch(choose_cort, nc);                                     %choose cortical patches of interest

%% Load Whitened Subcortical SVDs
ssvds_fname = strcat(fwd_path, prefix, '_subc_ico3_meg_', chtype, '_', gradmag, '_SVDs_rearranged.mat');
load(ssvds_fname,'sdiv_fwd','sregname');                                    %subcortical fwds by subdivisions with regnames
sregnums = 1:length(sregname);                                              %sregion indices to consider - excluded hipvols, else [2:5, 7:13]

%% Choose, Read and Merge Forward Solutions 
[Gtheta, Gn, currstr, lk_flag, lk, svolume, clkmax, reg_names] = cort_scort_fwds(lk, ico, cp, cort_cd, ...
sdiv_fwd, sregname, sregnums, nummodes, choose_scort);                      %store fwd structures for future reference

%% SRAS Coordinates and Pairwise Distances for all Sources
cfwd_fname = strcat(fwd_path, prefix,'-', meastype,'-all-fixed-fwd.fif');   %cortical forward solution
csrc_fname = strcat(prefix, '-ico-',num2str(cdiv),'p-src.fif');             %cortical source space
cfwd = mne_read_forward_solution(cfwd_fname);                               %cort fwd
csrc = mne_read_source_spaces(csrc_fname);                                  %cort source space
hfwd_fname = strcat(prefix, '_hipsurf-meg-all-fixed-fwd.fif');              %hippocampal forward solution
hsrc_fname = strcat(prefix, '_hipsurf-ico-1p-src.fif');                     %hippocampal source space
hfwd = mne_read_forward_solution(hfwd_fname);                               %hip fwd
hsrc = mne_read_source_spaces(hsrc_fname);                                  %hip source space
sras_coord = reg_coords(lk_flag, lk, clkmax, sdiv_fwd, hsrc, hfwd, csrc, cfwd); %sras coordinates
pw_dist = pairwise_dist(sras_coord);                                        %compute pairwise distances

%% Compute Whitened and Regularized Inverse Operator, and Resolution Matrix
C = eye(N,N);                                                               %whitened noise covariance matrix
sigmaX = N/trace(Gtheta*Gtheta');                                           %SNR of whitened data is 1/sigmaX
R0 = sigmaX*eye(size(Gtheta,2));                                            %source covariance matrix
lambda_sq = 1/SNR;                                                          %regularization parameter
W = (R0*Gtheta')/(Gtheta*R0*Gtheta'+lambda_sq*C);                           %inverse operator
K = W*Gtheta;                                                               %resolution matrix

%% Save and Plot
savepath = strcat(pwd, '/2016_Results_MatsFigs/');
savename = strcat(savepath, prefix, '_', choose_cort, '_cort_', choose_scort, '_scort_resmat.mat');
save(savename, 'prefix', 'datapath','choose_cort','choose_scort',...        %setup of resolution matrix problem
    'csvds_fname','cdiv', 'cfwd', 'ssvds_fname',...                         %where to find fwds and what params to use
    'Gtheta','reg_names','currstr','lk', 'clkmax','lk_flag',...             %forward solution, indices, current strengths
    'nummodes','SNR','sras_coord','pw_dist','C','R0','W','K');              %resolution matrix and summ. stats vars