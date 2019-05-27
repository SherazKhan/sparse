%% Data Paths
dpath = '/autofs/eris/purdongp/users/pavitrak/sourceloc_data/';
datapath = strcat(dpath, datatype, '/', prefix, '/');
addpath(genpath(datapath));
rmpath(genpath(strcat(datapath, 'fwd/'))); 
addpath(strcat(datapath, 'fwd/'));

%% Code Paths
addpath(genpath('/autofs/eris/purdongp/users/pavitrak/sourceloc_scripts/ScSPIGH'));
