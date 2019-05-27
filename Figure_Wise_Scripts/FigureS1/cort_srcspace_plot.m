%% Set Paths and Load Forward Solutions
clear; clc; close all;
prefix = 'SM04';                                                            %SM04, cat_004, manny
meastype = 'meg';                                                           %meg or meg_eeg
datatype = 'Illustrations';                                                 %AEP, SEP, Illustrations
datapath = strcat('/autofs/eris/purdongp/users/pavitrak/sourceloc_data/',datatype,'/',prefix);
g = 1; m = 25;                                                              %scaling ratios for gradiom vs. magnetom
gradmag = strcat(num2str(g),'g',num2str(m),'m');                            %naming string
chtype = 'all';                                                             %all channels

%% Load Angles
fname = strcat(datapath,'/fwd/', prefix,'_forangles_',meastype,'_',chtype,'_',gradmag,'.mat');
load(fname,'cpatch','lpatch_ind','sdiv','csvds_fname','ssvds_fname');
csvds_fname = strcat(datapath,'/fwd/', csvds_fname(find(csvds_fname=='/',1,'last')+1:end));
load(csvds_fname,'source','ico');
cfwd = mne_read_forward_solution(strcat(datapath, '/fwd/',prefix,'-', meastype,'-all-fixed-fwd.fif'));

for ico_num = 1:3;                                                          %ico division
    num_patch = length(ico(ico_num).patch);
    patchno = 1:1:num_patch; 
    plot_cort_surfacesources(source, ico, ico_num, patchno, cfwd, prefix);  %save a cortical patches surface map
end