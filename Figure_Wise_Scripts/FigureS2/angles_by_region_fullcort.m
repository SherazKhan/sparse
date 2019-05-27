%% Compute All Pairwise Angles
%addpath('/autofs/eris/purdongp/users/pavitrak/sourceloc_scripts/SM04-illust-code/angles_final');
%pairwise_angles_computation; pause; clc; close all;

%% Set Paths and Load Forward Solutions
clear; clc; close all;
prefix = 'SM04';                                                            %SM04, cat_004, manny
meastype = 'meg';                                                           %meg or meg_eeg
add_allpaths; datapath = strcat(datapath,'/');
ico_num = 2;                                                                %ico division (for visualization only)
chtype = 'all';                                                             %all channels
g = 1; m = 25;                                                              %scaling ratios for gradiom vs. magnetom
gradmag = strcat(num2str(g),'g',num2str(m),'m');                            %naming string

%% Load Angles
% use scort divisions sized to equal current strength of ico3 cort
fname = strcat(datapath,'fwd/', prefix,'_s2_forpairangles_',meastype,'_',chtype,'_',gradmag,'.mat');	
load(fname,'cpatch','sdiv','angfcs_ov','angcfs_ov','angcc_ov','angcs_ov','angss_ov','csvds_fname','ssvds_fname');
csvds_fname = strcat(datapath,'fwd/', csvds_fname(find(csvds_fname=='/',1,'last')+1:end));
ssvds_fname = strcat(datapath,'fwd/', ssvds_fname(find(ssvds_fname=='/',1,'last')+1:end));

angfcs = [];
for i= 1:length(angfcs_ov)
    angfcs = [angfcs angfcs_ov{i}'];
end
disp(max(angfcs))