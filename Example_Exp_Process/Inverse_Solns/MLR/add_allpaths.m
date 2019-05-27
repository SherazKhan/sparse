%% Data Paths
dpath = '/eris/purdongp/users/pavitrak/sourceloc_data';
switch datatype
    case 'AEP'
        datapath = strcat(dpath, '/', datatype, '/', prefix,'/');
        addpath(genpath(datapath));
        rmpath(genpath(strcat(datapath, 'meas')));
        addpath(strcat(datapath, 'meas/abr_plr_calibration'));
        addpath(strcat(datapath, 'meas/abr_plr_runall'));
    case 'SEP'
        datapath = strcat(dpath, '/', datatype, '/', prefix,'/');
        addpath(genpath(datapath));
    case 'Alpha'
        datapath = strcat(dpath, '/AEP/', prefix,'/');
        addpath(genpath(datapath));
end

%% Code Paths
codepath = strcat('/eris/purdongp/users/pavitrak/sourceloc_scripts/',prefix,'-code/');
addpath(codepath);
%addpath(genpath(strcat(dpath, '/', prefix,'inv-',prefix,'_',datatype)));
addpath(genpath('/eris/purdongp/users/pavitrak/sourceloc_scripts/ScSPIGH'));