clear; clc; close all; 
prefix = 'SM04';
sregions = {'lvpl','lstr','lhip','ramy','bsred'};                           %scenario name

for i = 1:5
    
srcname = sregions{i};
disp(srcname);
figs_path = strcat(pwd,'/random-',srcname,'_circ_results');
savepath = figs_path; 
saveon = 1;

load(strcat(figs_path,'/',prefix,'_random-',srcname,'_3c_angles_shuffle.mat'),'circ_name','angs_allcombos_sc','angs_allcombos_cc','num_cort');
figsavename = strcat(savepath, circ_name.overall,'_',num2str(num_cort), 'c_allmodecombos_shuffle');
[summ_stat_cc, summ_stat_sc] = figure_rdraws_results(angs_allcombos_sc, angs_allcombos_cc, saveon, figsavename);
disp('3c_shuffle')
size(angs_allcombos_sc)
size(angs_allcombos_cc)
clear circ_name angs_allcombos_sc angs_allcombos_cc num_cort;

load(strcat(figs_path,'/',prefix,'_random-',srcname,'_4c_angles_shuffle.mat'),'circ_name','angs_allcombos_sc','angs_allcombos_cc','num_cort');
figsavename = strcat(savepath, circ_name.overall,'_',num2str(num_cort), 'c_allmodecombos_shuffle');
[summ_stat_cc, summ_stat_sc] = figure_rdraws_results(angs_allcombos_sc, angs_allcombos_cc, saveon, figsavename);
disp('4c_shuffle')
size(angs_allcombos_sc)
size(angs_allcombos_cc)
clear circ_name angs_allcombos_sc angs_allcombos_cc num_cort;

load(strcat(figs_path,'/',prefix,'_random-',srcname,'_5c_angles_shuffle.mat'),'circ_name','angs_allcombos_sc','angs_allcombos_cc','num_cort');
figsavename = strcat(savepath, circ_name.overall,'_',num2str(num_cort), 'c_allmodecombos_shuffle');
[summ_stat_cc, summ_stat_sc] = figure_rdraws_results(angs_allcombos_sc, angs_allcombos_cc, saveon, figsavename);
disp('5c_shuffle')
size(angs_allcombos_sc)
size(angs_allcombos_cc)
clear circ_name angs_allcombos_sc angs_allcombos_cc num_cort;

end