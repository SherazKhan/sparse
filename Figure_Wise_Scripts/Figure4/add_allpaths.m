codepath = strcat(dpath, 'sourceloc_scripts/');
datapath = strcat(dpath, 'sourceloc_data/', datatype,'/',prefix);
addpath(genpath(strcat(codepath, 'ScSPIGH/')));
addpath(strcat(datapath, '/fwd/'));
addpath(strcat(datapath, '/src/')); 
addpath(strcat(datapath, '/mri/'));
addpath(strcat(datapath, '/meas/'));