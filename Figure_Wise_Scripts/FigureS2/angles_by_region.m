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
load(fname,'cpatch','lpatch_ind','sdiv','angfcs_ov','angcfs_ov','angcc_ov','angcs_ov','angss_ov','csvds_fname','ssvds_fname');
csvds_fname = strcat(datapath,'fwd/', csvds_fname(find(csvds_fname=='/',1,'last')+1:end));
ssvds_fname = strcat(datapath,'fwd/', ssvds_fname(find(ssvds_fname=='/',1,'last')+1:end));

%% Likely Angles by Region
% use scort divisions sized to equal current strength of ico-3 cort but do angles with ico-2 cort for visualization efficacy (no diff otherwise)
savename = strcat('angles_by_region_','_',chtype,'.mat');                   %Save Likely Angles in this File
save_plots = 1; ang_flag = 'median';        chk = 0;                        %Type of summary statistic
%save_plots = 1; ang_flag = 'mode';         chk = 1;%
if exist(savename,'file') ~=2
    disp('each cort vs. every other cort');     mkdir cort_vs_cort_hists;         
    lik_ccang = plot_cc_byregion(angcc_ov, cpatch, save_plots, ang_flag);          pause; close all;
    disp('each scort vs. every cort');          mkdir scort_vs_cort_hists;   
    slik_csang = plot_cs_byregion_s(angcs_ov, sdiv, cpatch, save_plots, ang_flag); pause; close all;
    disp('each cort vs. every scort');          mkdir cort_vs_scort_hists;     
    clik_csang = plot_cs_byregion_c(angcs_ov, sdiv, cpatch, save_plots, ang_flag); pause; close all;
    disp('each scort vs. every other scort');   mkdir scort_vs_scort_hists;        
    lik_ssang = plot_ss_byregion(angss_ov, sdiv, save_plots, ang_flag);            pause; close all;
    save(savename,'cpatch','lpatch_ind','sdiv','lik_ccang','slik_csang','clik_csang','lik_ssang');
end
load(savename);
if chk == 0                                                                 %Median Case
    lik_ccang(isnan(lik_ccang)) = 0;
    slik_csang(isnan(slik_csang)) = 0;
    clik_csang(isnan(clik_csang)) = 0;
    lik_ssang(isnan(lik_ssang)) = 0;
elseif chk == 1                                                             %Mode Case
    lik_ccang(lik_ccang == chk) = 0; 
    slik_csang(slik_csang == chk) = 0;
    clik_csang(clik_csang == chk) = 0;
    lik_ssang(lik_ssang == chk) = 0;
end

%% Plot Cortical Topographs
load(csvds_fname,'source','ico');
num_patch = length(ico(ico_num).patch);
patchno = 1:1:num_patch; 
cfwd = mne_read_forward_solution(strcat(prefix,'-', meastype,'-all-fixed-fwd.fif'));
plot_cort_angletopo(source, ico, ico_num, patchno, cfwd, prefix, [], lik_ccang, strcat('cc_',chtype));
plot_cort_angletopo(source, ico, ico_num, patchno, cfwd, prefix, [], clik_csang, strcat('cs_',chtype));

%% Plot Subcortical Topographs
mgz_file_list = dir(strcat(datapath, '/mri/svolume_decomps/*.mgz'));
for i = 1:length(mgz_file_list)
	names_submask{i} = mgz_file_list(i).name;
end
hmgz_fname = names_submask(find(~cellfun(@isempty, strfind(names_submask,'hipsurf'))));
plot_scort_angletopo(prefix, datapath, sdiv, hmgz_fname, lik_ssang, strcat('ss_',chtype));
plot_scort_angletopo(prefix, datapath, sdiv, hmgz_fname, slik_csang, strcat('sc_',chtype));

%% MATLAB Plots of Distributions
close all; num_bins = 15;
figure, set(gcf,'color','white'); lik_ccang = lik_ccang(lik_ccang>0); 
[~, binlocs_1] = hist(lik_ccang,num_bins); 
d = diff(binlocs_1)/2;  edges_1 = [binlocs_1(1)-d(1), binlocs_1(1:end-1)+d, binlocs_1(end)+d(end)]; edges_1(2:end) = edges_1(2:end)+eps(edges_1(2:end));
histogram(lik_ccang, edges_1, 'Normalization', 'probability');
%[~,X] = hist(lik_ccang, num_bins); histnorm(lik_ccang, X); clear X; 
xlim([0 90]); %axis([0 90 0 0.2]);
title('Cort vs. Cort'); ylabel('Angle (Degrees)'); xlabel('Region Index');
print(gcf,'-depsc','CC_Likely_Angles');

figure, set(gcf,'color','white'); lik_ssang = lik_ssang(lik_ssang>0);
[~, binlocs_2] = hist(lik_ssang,num_bins); 
d = diff(binlocs_2)/2;  edges_2 = [binlocs_2(1)-d(1), binlocs_2(1:end-1)+d, binlocs_2(end)+d(end)]; edges_2(2:end) = edges_2(2:end)+eps(edges_2(2:end));
histogram(lik_ssang, edges_2, 'Normalization', 'probability');
%[~,X] = hist(lik_ssang, num_bins); histnorm(lik_ssang, X); clear X; 
xlim([0 90]); %axis([0 90 0 0.2]);
title('Scort vs. Scort'); ylabel('Angle (Degrees)'); xlabel('Region Index');
print(gcf,'-depsc','SS_Likely_Angles');

figure, set(gcf,'color','white'); clik_csang = clik_csang(clik_csang>0);
[~, binlocs_3] = hist(clik_csang,num_bins); 
d = diff(binlocs_3)/2;  edges_3 = [binlocs_3(1)-d(1), binlocs_3(1:end-1)+d, binlocs_3(end)+d(end)]; edges_3(2:end) = edges_3(2:end)+eps(edges_3(2:end));
histogram(clik_csang, edges_3, 'Normalization', 'probability');
%[~,X] = hist(clik_csang, num_bins); histnorm(clik_csang, X); clear X; 
xlim([0 90]);% axis([0 90 0 0.2]);
title('Cort vs. Scort'); ylabel('Angle (Degrees)'); xlabel('Region Index');
print(gcf,'-depsc','CS_Likely_Angles');

figure, set(gcf,'color','white'); slik_csang = slik_csang(slik_csang>0);
[~, binlocs_4] = hist(slik_csang,num_bins); 
d = diff(binlocs_4)/2;  edges_4 = [binlocs_4(1)-d(1), binlocs_4(1:end-1)+d, binlocs_4(end)+d(end)]; edges_4(2:end) = edges_4(2:end)+eps(edges_4(2:end));
histogram(slik_csang, edges_4, 'Normalization', 'probability');
%[~,X] = hist(slik_csang, num_bins); histnorm(slik_csang, X); clear X; 
xlim([0 90]); %axis([0 90 0 0.2]);
title('Scort vs. Cort'); ylabel('Angle (Degrees)'); xlabel('Region Index');
print(gcf,'-depsc','SC_Likely_Angles');