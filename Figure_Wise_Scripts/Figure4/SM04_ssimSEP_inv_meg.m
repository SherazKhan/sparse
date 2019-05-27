%% Initialization
clear; clc; close all;
prefix = 'SM04';                                                            % Subject: SM04, cat004, manny
meastype = 'meg';                                                           % Measurement Type: meg or meg-eeg
datatype = 'Illustrations';                                                 % Data Paradigm
dpath = '/autofs/eris/purdongp/users/pavitrak/';  add_allpaths;             % Add Paths
g = 1; m = 25; gradmag = strcat(num2str(g), 'g',num2str(m),'m'); chtype = 'all'; % grad mag ratio
%diary(strcat(prefix,'_',meastype,'_matlog'));                              % log the run
savepath = strcat(pwd,'/2016_Results_MatsFigs/'); 

%% Settings for Problem Setup and Greedy Estimation
settings_fname = strcat(datapath, '/fwd/',prefix,'_ico3_forangles_meg_all_',gradmag,'.mat');
load(settings_fname,'Cw','meg_sel'); sel = meg_sel; N = length(sel); clear meg_sel;% Read whitener (diagonal, preserves order) & MEG chs
cij = zeros(1,N); cij(1:3:N) = 1;                                       % Flag Gradiometers vs. Magnetometers 
alg = 'MNE'; stop = [];                                                     % Parameters for Greedy

%% Sparse Cortical Patch Info and Subcortical Info
load(strcat(savepath, prefix, '_sparse_cort_all_scort_resmat.mat'),...
    'csvds_fname','cdiv','lk','clkmax', 'ssvds_fname','nummodes','SNR');
patch_est = lk(lk <= clkmax);                                               %selected cortical patches

%% Save Important Variables and Call Inverse Solutions
savename = strcat(savepath, prefix,'-',meastype,'-',datatype,'_sinv.mat');
save(savename, 'prefix','meastype','datatype','csvds_fname','ssvds_fname',...
               'Cw','sel','alg','stop','cij','nummodes','SNR','cdiv','patch_est','-v7.3');
SM04_gen_noiseless_sim;
SM04_final_scort_greedy_meg;
%diary off;