%% Set Paths and Load Forward Solutions
%General Paths and Subject Information
clear; clc; close all;
addpath(genpath('/autofs/eris/purdongp/users/pavitrak/sourceloc_scripts/ScSPIGH'));
prefix = 'SM04'; meastype = 'meg'; datatype = 'Illustrations';
subjdir = strcat('/autofs/eris/purdongp/users/pavitrak/sourceloc_data/',datatype, '/',prefix);
chtype = 'all';                                                             % Consider both magnetometers and gradiometers
fwd_path = strcat(subjdir, '/fwd/');                                        % Forward solutions are all stored here
mri_path = strcat(subjdir,'/mri/');                                         % MRI path for display of spatial distribution of current estimates
g = 1; m = 25;                                                              % Scaling ratios for gradiom vs. magnetom
gradmag = strcat(num2str(g),'g',num2str(m),'m');                            % Naming string

% USER SPECIFY Key Parameters                                               %%%%%%%%%%%%%%%%% USER INPUTS %%%%%%%%%%%%%%%% 
scort_ind = 1;                                                              %1 for top mode, 7 for all modes for VPL
type_flag = 'all_top';                                                      %'all_top, 'all', 'top' 

%Field Map  Comparisons Parameters
circ_choice = 'som_erp'; disp(strcat('circuit # ',circ_choice));            %'inf_hip'; %specify circuit being analyzed
savepath = strcat(pwd,'/',circ_choice,'_circ_results/'); 
grp_fname = strcat(savepath, prefix,'_forgrpangles_',circ_choice,'_meg_all_',gradmag,'_sc.mat');
load(grp_fname, 'scortgrp','cortgrp');                                      %file that contains angles and subsets for comparisons
sdiv = scortgrp(1).sdiv; cpatch = cortgrp(1).cpatch;                        %subcortical and cortical structures
thresh = [30 45 60 80]; num_thresh = length(thresh);                        %threshold angles to pick mode subsets - 43.36, 43.5-45.2

%Fitting Parameters
ico_num = 2;                                                                % ICO-2 Cortical Patches
sregname = 'lt';  sregsdiv = 1;                                             % Left thalamus Div 1 for 2014 Thesis Version of VPL
[c_currstr, s_currstr] = currstr_calc(subjdir, prefix, meastype, chtype, gradmag, ico_num, sregname, sregsdiv); % get current strengths for VPL, all cortical patches
srctype = 'scort'; fittype = 'sparse_cort';                                 % Specify fit type for all output files saved
SNR_fit = 9;                                                                % Used to determine lambda_sq = 1/SNR - previously set to 1

%Field Map Plotting Parameters
load(strcat(fwd_path,prefix,'_ico2_forangles_meg_all_',gradmag,'.mat'),'Cw');    % Read in whitener (diagonal, preserves order) and MEG channels
phys_meas_fname = strcat(subjdir, '/meas/',prefix,'_SomSens_SingleNoTask_offl_SSP.fif'); % For plotting field maps
virt_meas_fname = strcat(subjdir, '/meas/',prefix,'_virt_meg.fif');         % For plotting field maps
helmet_fname = '/autofs/eris/purdongp/users/pavitrak/sourceloc_scripts/ScSPIGH/fieldmaps_angles/Elekta_VirtMag304.mat';
proj_approach = 'cart_sph';                                                 % Choose 'cart_sph' or 'mdimsc
num_contours = 0;                                                           % Choose # of contours to display - was 50

%% Identify Mode Subsets for Field Map Comparison
[~, angf1f2, sdiv_subsets, sfwd_emode, sfwd_eval, cpatch_subsets, cfwd_emode, cfwd_eval] ...
    = reg1_reg2_allcombos(sdiv, cpatch, type_flag, sdiv.sp, cpatch(1).cp);                         
num_angles_cort = size(angf1f2,1); num_angles_scort = size(angf1f2, 2);     %# cortical and subcortical options
st = 1:num_angles_cort:num_angles_scort*num_angles_cort;                    %starting point within angf1f2
en = num_angles_cort: num_angles_cort:num_angles_scort*num_angles_cort;     %ending point within angf1f2
ind = cell(1,num_thresh); subset_inds = struct([]);
for cnt = 1:num_thresh
  [~, cort_ind] = min(abs([angf1f2{st(scort_ind):en(scort_ind)}] - thresh(cnt)));
  subset_inds(cnt).thresh = thresh(cnt);
  subset_inds(cnt).cort = cort_ind;
  subset_inds(cnt).scort = scort_ind;
  subset_inds(cnt).angles = angf1f2{subset_inds(cnt).cort, subset_inds(cnt).scort};
  subset_inds(cnt).smodes = sdiv_subsets{subset_inds(cnt).scort};
  subset_inds(cnt).sfwd_emode = sfwd_emode(:,subset_inds(cnt).smodes);
  subset_inds(cnt).sfwd_eval = sfwd_eval(subset_inds(cnt).smodes); 
  subset_inds(cnt).cmodes = cpatch_subsets{subset_inds(cnt).cort};
  subset_inds(cnt).cfwd_emode = cfwd_emode(:,subset_inds(cnt).cmodes);
  subset_inds(cnt).cfwd_eval = cfwd_eval(subset_inds(cnt).cmodes);
end
save(strcat(savepath, 'fmap_settings.mat'),'angf1f2','subset_inds','srctype', 'fittype', 'scort_ind','gradmag', 'num_contours','proj_approach');

%% Initialize Comparison Variables
xs = cell(1,num_thresh);   Gs = xs; xfit = xs;     Gfit = xs;  ys = xs;    lkmax = zeros(1,num_thresh);
s_fmap_fname = xs;         c_fmap_fname = xs;      news_fmap_fname = xs;   newc_fmap_fname = xs; 
s_data_2D = xs;            s_data_2Dtoplot = xs;   c_data_2D = xs;         c_data_2Dtoplot = xs;

for i = 1:num_thresh
if ~isempty(subset_inds(i).smodes)
%% Compare Field Maps
[xs{i}, Gs{i}, xfit{i}, Gfit{i}, lkmax(i), ys{i}] = fit_scort_to_sparsecort(sdiv, subset_inds(i), SNR_fit, s_currstr, c_currstr);

%% For Field Maps: DeWhiten to Physical Units
Gs_no_whiten{i} = Cw^(-1)*Gs{i};                                            % Dewhitening Gs implicitly dewhitens ys used for field maps
Gfit_no_whiten{i} = Cw^(-1)*Gfit{i}(:,1:lkmax(i));                          % Dewhitening Gfit implicitly dewhitens ycfit used for field maps

%% Store For Field Maps (3D Physical Helmet)
% Field Map settings
data_3D = fiff_read_evoked_all(phys_meas_fname);                            % Read data file to modify
[chnames, chsel] = channel_selector('meg', data_3D, chtype);                % Select channels relevant to meastype
data_3D = fiff_pick_channels_evoked(data_3D, chnames);                      % Read data specifically for these channels
% Subcortical Activity Map
s_fmap_fname{i} = strcat(savepath, 'real_',strcat(srctype,num2str(i)),'_3D.fif');     % field map file name
fieldmap3D(data_3D, Gs_no_whiten{i}, s_fmap_fname{i}, chsel, xs{i});        % Write field map for mne_analyze view      
% Cortical Fit Map
c_fmap_fname{i} = strcat(savepath,'real_',strcat(fittype,num2str(i)),'_3D.fif'); % field map file name
fieldmap3D(data_3D, Gfit_no_whiten{i}, c_fmap_fname{i}, chsel, xfit{i}(1:lkmax(i),:)); % Write field map for mne_analyze view      

%% Convert Field Maps from Physical Helmet to Virtual Magnetometer Array 
news_fmap_fname{i} = strcat(savepath, 'virt_',strcat(srctype,num2str(i)),'_3D.fif');  % field map file name
newc_fmap_fname{i} = strcat(savepath, 'virt_',strcat(fittype,num2str(i)),'_3D.fif');  % field map file name
disp('run the following in terminal for forward solutions');
clear unix_cmd;
unix_cmd{1} = ['mne_map_data --from ', s_fmap_fname{i}, ' --frommri $SUBJECTS_DIR/$SUBJECT/fwd/$SUBJECT"-trans.fif" --to ',...
                virt_meas_fname, ' --tomri $SUBJECTS_DIR/$SUBJECT/fwd/$SUBJECT"-trans.fif" --out ', news_fmap_fname{i}];
unix_cmd{2} = ['mne_map_data --from ', c_fmap_fname{i}, ' --frommri $SUBJECTS_DIR/$SUBJECT/fwd/$SUBJECT"-trans.fif" --to ', ...
                virt_meas_fname, ' --tomri $SUBJECTS_DIR/$SUBJECT/fwd/$SUBJECT"-trans.fif" --out ', newc_fmap_fname{i}];
disp(char(unix_cmd(~cellfun(@isempty,unix_cmd))));
disp('after it completes: type exit at terminal');                          % Should be Fast
pause; clc;

%% Plot Field Maps (2D Virtual Magnetometer)
% Subcortical Activity Map
s_data_2D{i} = fiff_read_evoked_all(news_fmap_fname{i});                    % subcortical data in virtual mag system
s_data_2Dtoplot{i} = s_data_2D{i}.evoked.epochs(:,1);                       % only 1 time point required
s_scale(i) = max(abs(min(s_data_2Dtoplot{i})), abs(max(s_data_2Dtoplot{i})));
s_data_2Dtoplot_sc{i} = s_data_2Dtoplot{i}/s_scale(i);
s_ax = fieldmap2D(helmet_fname, proj_approach, virt_meas_fname, s_data_2Dtoplot_sc{i}, num_contours,[]);
set(gcf,'renderer','painters'); print(gcf,'-depsc',strcat(savepath, 'fmap_', strcat(srctype,num2str(i)),'_src'));% -depsc for color
% Cortical Fit Map
c_data_2D{i} = fiff_read_evoked_all(newc_fmap_fname{i});                    % cortical fit in virtual mag system
c_data_2Dtoplot{i} = c_data_2D{i}.evoked.epochs(:,1);                       % only 1 time point required
c_scale(i) = max(abs(min(c_data_2Dtoplot{i})), abs(max(c_data_2Dtoplot{i})));
c_data_2Dtoplot_sc{i} = c_data_2Dtoplot{i}/c_scale(i);
c_ax = fieldmap2D(helmet_fname, proj_approach, virt_meas_fname, c_data_2Dtoplot_sc{i}, num_contours,s_ax);
set(gcf,'renderer','painters'); print(gcf,'-depsc',strcat(savepath, 'fmap_',strcat(fittype,num2str(i)),'_fit')); % -depsc for color
clear s_ax c_ax;

%% Compute Metrics to Report
numer = norm(s_data_2Dtoplot{i} - c_data_2Dtoplot{i}).^2;                   % numerator for GOF  - MSE
denom = norm(s_data_2Dtoplot{i}).^2;                                        % denominator for GOF - MSE
gof = 1 - sqrt(numer/denom);                                                % GOF b/w Subcortical Source & Best Fit Cortical Field Maps
pang = max(subspacea(Gs{i}, Gfit{i}))*180/pi;                               % principal angle for this specific combination
title({['angle  ', num2str(subset_inds(i).angles)], ['cortical index  ', num2str(subset_inds(i).cort)], ...
       ['max ang    ', num2str(pang)],  ['   GOF  ', num2str(gof)]});       % cortical index
end
end
disp(c_scale);