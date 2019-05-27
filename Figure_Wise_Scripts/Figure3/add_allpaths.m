%% Data Paths
dpath = '/eris/purdongp/users/pavitrak/sourceloc_data/Illustrations';
datapath = strcat(dpath, '/',prefix,'/');
addpath(genpath(datapath));
rmpath(strcat(datapath, 'fwd/2015_VirtualMag_System'));
rmpath(strcat(datapath, 'fwd/2014-15_Vectorview_System'));

%% Code Paths
addpath(genpath('/eris/purdongp/users/pavitrak/sourceloc_scripts/ScSPIGH'));