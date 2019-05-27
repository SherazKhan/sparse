%%%%%%%%%%%%%%%%%%Analysis of Resolution Matrix and Plots%%%%%%%%%%%%%%%%%%
% Written by Pavitra Krishnaswamy (pavitrak@nmr.mgh.harvard.edu
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Paths and Parameters
clear; clc; close all;
prefix = 'SM04'; meastype = 'meg'; datatype = 'Illustrations';
dpath = '/autofs/eris/purdongp/users/pavitrak/'; add_allpaths;               %add all paths
datapath = strcat(dpath, 'sourceloc_data/', datatype,'/',prefix,'/');
add_allpaths;
saveon = 1; savepath = strcat(pwd, '/2016_Results_MatsFigs/');
flag_on = 1;

for testcase = 1:3
%% Assign Files Containing Resolution Matrix and Centroid Distance Information
[resmat_name, cax, dist_name, savefigname, savestatsname] = assign_type(testcase, prefix, savepath);
load('revised_Fig4_Colormap'); cmap  = eval(strcat('cmap',num2str(testcase)));

%% Plot Resolution Matrix
load(resmat_name,'K','nummodes','SNR'); load(dist_name,'reg_names','cdiv','lk','clkmax');
disp(['nummodes:  ', num2str(nummodes)]);
disp(['SNR for regularization:  ', num2str(SNR)]); 
reg_list = unique(reg_names,'stable');
for rr = 1:length(reg_list)
    st(rr) = find(strcmp(reg_names, reg_list{rr}), 1, 'first');
    en(rr) = find(strcmp(reg_names, reg_list{rr}), 1, 'last');
end
reg_matrix_plots(abs(K), reg_list, st, en, saveon, savefigname, cax, flag_on, cmap); shg; %pause

%% Pairwise Distances
if testcase == 3  %&& nummodes == 0                                          % Greedy Case
    load(dist_name, 'sras_coord','lk_flag');
    new_lk_flag = ones(1,sum(lk<=clkmax)); clk = lk(lk<=clkmax);
    cfwd = mne_read_forward_solution(strcat(prefix,'-', meastype,'-all-fixed-fwd.fif'));      
    csrc = mne_read_source_spaces(strcat(prefix, '-ico-',num2str(cdiv),'p-src.fif'));
    cort_coord = reg_coords(new_lk_flag, clk, [], [], [], [], csrc, cfwd);  % SRAS coordinates
    new_sras_coord = [cort_coord sras_coord(~lk_flag)];
    pw_dist = pairwise_dist(new_sras_coord);                                % Compute pairwise distances
    A = unique(pw_dist, 'rows','stable'); B = unique(A','rows','stable'); pw_dist = B';
else
    load(dist_name,'pw_dist');                                              % In Centimeter
end
    
%% Spatial Dispersion
for i = 1:size(K,2)
	sd_num(i,:) = (pw_dist(:,i).*K(:,i)).^2;
    sd_den(i,:) = (K(:,i)).^2;
    sd_num_sum(i) = sum(sd_num(i,:),2);
    sd_den_sum(i) = sum(sd_den(i,:),2);
    spatial_disp(i) = sqrt(sd_num_sum(i)./sd_den_sum(i));                   % In Centimeter
end

%% Dipole Localization Error
for i = 1:size(K,2)
    [~,jj_est] = max(abs(K(:,i)));                                          % Max. point of column i -- where peak est is
    dle(i) = pw_dist(jj_est, i);                                            % Bias in Centimeter
end

%% Plot Histograms of Summary Statistics 
figure, set(gcf,'color','white'); 
subplot(2,1,1), histogram(spatial_disp, 20, 'Normalization','probability'); 
xlabel('Spatial Dispersion (cm)'); ylabel('Histogram'); xlim([0 7]);
subplot(2,1,2), histogram(dle, 20, 'Normalization','probability');
xlabel('Dipole Localization Error (cm)'); ylabel('Histogram'); xlim([0 7]);
print('-depsc', savestatsname);
print('-dpdf',savestatsname);
   
figure, set(gcf,'color','white'); 
c_ind = lk<=clkmax;
subplot(2,1,1), hold all;
histogram(spatial_disp(~c_ind), 20, 'FaceColor','b','FaceAlpha',0.5);
histogram(spatial_disp(c_ind),20,'FaceColor','r','FaceAlpha',0.5); 
xlabel('Spatial Dispersion (cm)'); ylabel('Histogram'); xlim([0 7]); legend('Cortical','Subcortical');
subplot(2,1,2), hold all;
histogram(dle(~c_ind), 20, 'FaceColor','b','FaceAlpha',0.5);
histogram(dle(c_ind),20,'FaceColor','r','FaceAlpha',0.5); 
xlabel('Dipole Localization Error (cm)'); ylabel('Histogram'); xlim([0 7]); legend('Cortical','Subcortical');

%% Save Plots and Results
if saveon  
   savename = strcat(savepath, prefix, '_testcase',num2str(testcase), '_plot_summary.mat');
   save(savename, 'K', 'nummodes','SNR','reg_list', 'st', 'en', 'savefigname', 'cax', 'spatial_disp', 'dle', 'savestatsname'); 
end
clear resmat_name dist_name cax savefigname savestatsname K reg_names cdiv reg_list st en...
      pw_dist sd_num sd_den sd_num_sum sd_den_sum spatial_disp jj_est dle;
pause; %close all;
end