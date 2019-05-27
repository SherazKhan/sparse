%% Initialize Paths and Filenames
clear; clc; close all;
prefix = 'SM04';                                                            % SM04, cat_004, manny
meastype = 'meg';                                                           % meg or meg_eeg
datatype = 'SEP';                                                           % Illustrations, SEP, AEP
append = '_overlap';                                                        % study results for case with or w/o overlap
fullcort = 1;                                                               % 0 with hierarchy, 1 without hierarchy
cdiv = 3;  hfwd_type = 'ico-1';                                             % cortical and hippocampal source spaces
if strcmp(append,'overlap')
masks_tstep = 25;                                                           % time step for masks in milliseconds
greedy_mask_st = 9;                                                         % starting point in milliseconds
mne_mask_st = 59;                                                           % starting point in milliseconds
mne_mask_en = 109;                                                          % end point in milliseconds
else
masks_tstep = 9;                                                            % time step for masks in milliseconds
greedy_mask_st = 0;                                                         % starting point in milliseconds
mne_mask_st = 0;                                                            % starting point in milliseconds
mne_mask_en = 54;                                                           % end point in milliseconds
end
if fullcort 
    thresh = 90; thresh2 = 90; 
else
    thresh = 95; thresh2 = 100; 
end
time_norm = [];                                                             % compute norm for all times before time_norm (msec)

%% Include Paths, Specify File Names
add_allpaths;                                                               % add paths
if fullcort == 0
    scort_resfname = strcat(pwd,'/',append(2:end),'_final/', prefix,'-',meastype,'-',datatype,'_scortest',append,'.mat');       % subcortical results file
else
    scort_resfname = strcat(pwd,'/',append(2:end),'_final/', prefix,'-',meastype,'-',datatype,'_1stage_scortest',append,'.mat');% subcortical results file
end
meas_fname = strcat(datapath,'meas/', prefix,'_vpl_somato_meg_sim',append,'.mat'); % simulated measurement file
subc_svds_fname = strcat(datapath, 'fwd/', prefix,'_subc_',meastype,'_prewhite','_SVDs',append,'.mat'); %subcortical SVD file
cort_svds_fname = strcat(datapath, 'fwd/', prefix,'_cort-meg-all_white_SVDs',append,'.mat'); %cortical SVD file
curr_den_fname = strcat(datapath, 'fwd/curr_den.mat');                      %for cortical current density
diary(strcat(prefix, '_hybrid-meg-results_visuals'));

%% Load Measurement, Subcortical Params & Results, Fwd & Src, ASEG & T3
load(meas_fname,'Xc_plot','Xvpl','Xvpl_plot','scalarX','stc','t','sim_on'); alg = 'MNE'; % load inputs into subcortical inverse problem
load(scort_resfname);                                                       % load hybrid estimation results
cfwd = mne_read_forward_solution(strcat(datapath,'fwd/', prefix, '-meg-all-fixed-fwd.fif'));         %cort fwd
csrc = mne_read_source_spaces(strcat(datapath,'src/', prefix, '-ico-',num2str(cdiv),'p-src.fif'));   % cort source space
hfwd = mne_read_forward_solution(strcat(datapath,'fwd/', prefix, '_hipsurf-meg-all-fixed-fwd.fif')); % hip fwd
hsrc = mne_read_source_spaces(strcat(datapath,'src/', prefix, '_hipsurf-',hfwd_type,'p-src.fif'));   % hip source space
aseg = MRIread(strcat(datapath,'mri/aseg.mgz'));                            % aseg for mask creation
load(strcat(datapath,'/mri/T3.mat'));                                       % sRAS to vRAS transform
load(subc_svds_fname,'sdiv_fwd','sregname');                                %subcortical fwds by subdivisions with regnames
load(cort_svds_fname); eigmodes = p; clear p;                               %cortical fwds patches   

%% Identify Relevant Cortical vs. Subcortical Elements of Solution
% Indices for All Regions, Greedy Choices, and their Cortical and Subcortical Components
allind = 1:1:size(Gtheta,2);                                                % all indices in fwd solution
cindmne = lk <= clkmax; cmne = allind.*cindmne;                             % indicators for all cortical vs subcortical
sindmne = 1-cindmne; smne = allind.*sindmne;                                % divide all indices into cortical and subcortical pieces
cind = lk(trimP) <= clkmax; ctrimP = trimP.*cind;                           % indicators for selected cortical vs subcortical
suind = 1-cind; strimP = trimP.*suind;                                      % divide trimP into cortical and subcortical pieces
for i = 1:length(regname)
    ind_sfwd{i} = strncmp(regname{i}, all_regions, 4);                      % sdiv_fwd indices for ith scort region in srcspace            
end
[hipsurf_ind1, reg_index, curr_str, curr_str_greedy, cpatch, sdiv_fwd_plot, all_regions] = inds2modes(sdiv_fwd, regname, ...
         all_regions, lk, trimP, ind_sfwd, ssubdiv, svolume, currstr, snummodes, clkmax, cp); % Index all variables by # Modes
clear suind count pw ind_sfwd sdiv_fwd sdiv_srcpspace0 sdiv_fwd_plot0 sdiv_srcspace ind_hipsurfs ...
      C dub N select_time allchannel sregnums subcort_ind t1 t2 i 

%% Read Measurement, Estimates (Already Scaled by Current Strengths)
Xls_scaled = Xgreedy{end};                                                  % greedy source current estimates for modes
Xls_mne_scaled = Xmne;                                                      % prescaled mne source current estimates for modes
T = size(Xls_scaled,2);                                                     % length of time period under study
Xsim = Xvpl(1:T); Xsim_plot = Xvpl_plot;                                    % simulated source current
time_est = time_est*1000;                                                   % milliseconds
clear Xgreedy Xmne;                                                         % make room for memory intensive tasks

%% Identify Regions and Assign Labels for Modes Found
[leg_fig, clabel_fname, hlabel_fname] = ID_legends_labels(lk, trimP, clkmax, strimP, ctrimP, csrc, cfwd, ...
hsrc, hfwd, hipsurf_ind1, regname, ssubdiv, svolume, datapath, prefix, cdiv); % ID Regions & Labels
pause; 
orig_rois = {'cS1','cS2','PPC','iS2','VPL'};                                % regions of interest
disp('please load labels or submasks into mne-analyze or freeview and confirm on corresponding roi')
disp(orig_rois);
if fullcort == 0
    if strcmp(append,'_overlap')
    rois = {'cs1','cs1','cs2','cs1','is2','cs1','ppc','vpl','cs2','rt','cs1','ppc'};%for overlap SM04 results
    else
    rois = {'cs1', 'cs1', 'ppc', 'is2', 'cs2', 'cs2', 'rt','cs2', 'is2', 'cs1', 'ppc','lc','vpl'}; %for no-overlap SM04 results
    end
else
    if strcmp(append,'_overlap')
    rois = {'cs2', 'cs1', 'cs1', 'other', 'ppc', 'cs1', 'cs2', 'is2', 'cs1', 'cs2', 'ppc', 'cs1'};
    else
    rois = {'cs2', 'cs2', 'cs1', 'is2', 'cs1', 'ppc', 'other1', 'cs1', 'cs1', 'other2', 'cs2', 'ppc', 'other3'};
    end
end
cort_or_scort = cind; disp(rois);
cortind = find(~cellfun(@isempty,clabel_fname));                            % cortical surface labels
[~,indc, ~] = unique(clabel_fname(cortind)); unique_patches = cortind(indc);% find the unique patch indices
scortind = find(cellfun(@isempty,clabel_fname));                            % subcortical volumes and hipsurfs
[~,inds, ~] = unique(leg_fig(scortind)); unique_vols = scortind(inds);      % find the unique volume indices
uni_divs = [unique_patches unique_vols];                                    % unique division indices in trimP
rois_gr = rois(uni_divs);                                                   % dont include rois for repeat patches or divs
cort_scort_gr = cort_or_scort(uni_divs);                                    % cortical or subcortical flag
[unique_rois, uni_inds, ~] = unique(rois_gr);                               % only consider unique region labels
unique_cort_scort = cort_scort_gr(uni_inds);                                % cortical or subcortical flag
col = hsv(length(unique_rois));                                             % colors for each ROI

%% Greedy & MNE BackProject and Plot MATLAB Time Courses
[Xls_plot_greedy, Xls_grind_mode] = backproj_greedy(lk, clkmax, trimP, strimP, ctrimP, reg_index, svolume, cpatch, sdiv_fwd_plot, ...
    ico(cdiv), Xls_scaled, curr_str, hsrc, hfwd, csrc); 
Xls_plot_grcond = roi_sum(rois_gr, unique_rois, uni_divs, Xls_plot_greedy,'greedy'); % condensed greedy time course estimates by roi
if fullcort
  for i = 1:size(Xls_mne_scaled, 1) 
      n(i) = norm(Xls_mne_scaled(i,:));
  end
  ind_sig_mne = find(n > prctile(n,thresh));
else
  ind_sig_mne = 1:length(allind);
end
[Xls_plot_mne_red, ~] = backproj_mne(lk, clkmax, allind(ind_sig_mne), smne(ind_sig_mne), cmne(ind_sig_mne), ...
    reg_index, svolume, cpatch, sdiv_fwd_plot, ico(cdiv), Xls_mne_scaled, curr_str, hsrc, hfwd, csrc); 
Xls_plot_mne(1:length(allind)) = struct('numdipoles',[],'cstrdip',[],'orient',[],'dip_or',[],'dip_res',[],...
        'dip_resor',[],'dip_resovern',[],'dip_resnorm',[], 'mode_indices',[],'mode',[]);    % initialize
for i = 1:length(allind)
    Xls_plot_mne(i).dip_res = zeros(1, T);                                  %initialize
end
Xls_plot_mne(ind_sig_mne) = Xls_plot_mne_red;                               % assign to variable we will use for plotting masks
Xls_plot_mnecond_gr = roi_sum(rois_gr, unique_rois, trimP(uni_divs), Xls_plot_mne,'mne-gr');  % same greedy rois to condense mne time course, to compare

%% MRI Movie Visualizations
v2ras1 = aseg.vox2ras1; mult = inv(v2ras1)*T3;                              % transform between sRAS and voxels

% Greedy & MNE MRI Masks/Voxel Coordinates - transform sRAS to voxels
disp('making greedy mask');
[mask_greedy, mask_offset_greedy, mask_mean_greedy] = masks_greedy(lk, clkmax, trimP, strimP, ctrimP, sdiv_fwd_plot, ...
Xls_plot_greedy, mult, aseg, hlabel_fname, clabel_fname);
disp('making mne mask');
[mask_mne, mask_offset_mne, mask_mean_mne, vox_coord] = masks_mne(lk, clkmax, allind, smne, cmne, ...
    sdiv_fwd_plot, Xls_plot_mne, hsrc, hfwd, csrc, cfwd, mult, aseg);              

% Greedy & MNE MRI Mask with Time Movies
tstep = round(masks_tstep/mean(diff(time_est))); [~, Tst] = min(abs(time_est-greedy_mask_st));% parameters downsampled time series for MRI greedy mask
disp('making greedy mask: time movie')
write_mask_movie(trimP, mask_greedy, tstep, Tst, T, Xls_plot_greedy, mask_offset_greedy, mask_mean_greedy, ...
aseg, fullcort, prefix, meastype, scalarX);                                 % create greedy mask
disp('making mne mask: time movie');
[~, Tst] = min(abs(time_est-mne_mask_st)); [~, Ten] = min(abs(time_est - mne_mask_en)); % parameters downsampled time series for memory efficient MRI MNE mask
if fullcort
write_mask_movie_all(ind_sig_mne, vox_coord(ind_sig_mne), tstep, Tst, Ten, Xls_plot_mne_red, mask_offset_mne(ind_sig_mne), ...
mask_mean_mne(ind_sig_mne), aseg, fullcort, prefix, meastype, scalarX);     % create MNE mask
else
write_mask_movie_all(allind, vox_coord, tstep, Tst, Ten, Xls_plot_mne, mask_offset_mne, mask_mean_mne, ...
aseg, fullcort, prefix, meastype, scalarX);                                 % create MNE mask
end

% Greedy & MNE Cortical Inflated Surfaces with Time Movies
ctrimP_nz = ctrimP(ctrimP>0); p = eigmodes(cdiv); 
patchno_gr = unique(lk(ctrimP_nz));                                         % patch nos. within ico-3 sourcespace
trimP_gr = unique(ceil(sort(ctrimP_nz)/p));                                 % index nos. within hybrid sourcespace
Xgr_hat = Xls_scaled(:,1:T);                                                % currents within hybrid sourcespace
Xc_plot_greedy = full_recon(alg,ico,cfwd,p,patchno_gr,prefix,source,stc,cdiv,trimP_gr,Xgr_hat,L,curr_den_fname,scalarX);
cmne_nz = cmne(ind_sig_mne); cmne_nz = cmne_nz(cmne_nz > 0);                % top MNE modes
patchno_mne = unique(lk(cmne_nz));                                          % patch nos. within ico-3 sourcespace
trimP_mne = unique(ceil(sort(cmne_nz/p)));                                  % index nos. within hybrid sourcespace
Xmne_hat = Xls_mne_scaled(:,1:T);                                           % currents within hybrid sourcepsace
Xc_plot_mne = full_recon(alg,ico,cfwd,p,patchno_mne,prefix,source,stc,cdiv,trimP_mne,Xmne_hat,clkmax,curr_den_fname,scalarX);

%% Plot Time Courses, Current Distributions Across Regions, Angles
% Sort out ROIS, Legends, Currents Strorage for MNE
if fullcort 
   cort_roi_index = num2cell(patchno_mne); 
   for i = 1:length(cort_roi_index)
       [clblno(i), ~, ~, clr(i)] = patch_to_labloc(csrc,cfwd,cort_roi_index{i});% label # corresponding to patch
       cort_rois{i} = strcat(num2str(clblno(i),'%06g'),'-',clr(i),'h.label');  %cortical label file name
   end
else
    if strcmp(append,'_overlap')
       cort_rois = {'ppc','ppc','cs1','cs1','cs1','cs2','cs2','is2'};
       cort_roi_index = {1,2,1,2,3,1,2,1};
   else
       cort_rois = {'ppc','cs1','cs1','cs2','cs2','cs2','is2'};
       cort_roi_index = {1, 1, 2, 1, 2, 3, 1};
   end
end
l = length(cort_rois); 
[uni_reg_index, uniq_regs, ~] = unique(reg_index(ind_sig_mne,:),'rows','stable');
legs_mne = cell(1, size(uni_reg_index,1)); flag_mne = zeros(1,size(uni_reg_index,1));
for i = 1:length(legs_mne)
    aa_chk = uni_reg_index(i,1) > length(regname); bb_chk = (fullcort == 1) && (i < length(cort_rois));
    if  aa_chk || bb_chk
        legs_mne_cnonum{i} = cort_rois{i};
        legs_mne{i} = strcat(cort_rois{i}, '-',num2str(cort_roi_index{i}));
        flag_mne(i) = 1;
    else
        legs_mne_snonum{i} = regname{uni_reg_index(i,1)};
        legs_mne{i} = strcat(legs_mne_snonum{i}, '-', num2str(uni_reg_index(i,2)));
    end
end
Xls_plot_mne_unique = Xls_plot_mne(ind_sig_mne(uniq_regs));                 % get rid of repeating mne time courses
l = min(l, length(Xls_plot_mne_unique));                                    % get rid of repeating mne time courses
if fullcort == 0
    [Xls_plot_cmnecond, legs_mne_crevise] = roi_sum(legs_mne_cnonum, unique(legs_mne_cnonum), allind(cindmne), ...
    Xls_plot_mne_unique(1:l), 'mne-cort', legs_mne(1:l));                   % condense mne time courses and legends for cortical
else
    Xls_plot_cmnecond = Xls_plot_mne_unique;
    legs_mne_crevise = legs_mne_cnonum;
    for i = 1:length(clabel_fname)
	ss = strcmp(cort_rois, clabel_fname(i)); sum(ss)
	if sum(ss) > 0
        ind(i) = find(ss);                                                  %find the roi indices
        cort_rois(ind(i)) = rois(i);                                        %replace cortical label
        legs_mne_crevise(ind(i)) = rois(i);                                 %replace cortical legends
	end
    end
end
Xls_plot_mnecond_all = [Xls_plot_cmnecond Xls_plot_mne_unique(l+1:end)];    %condensed mne time courses for all
legs_mne_condall = [unique(legs_mne_crevise) legs_mne(l+1:end)];            % legends for condensed time courses for all
cort_scort_mne = [ones(1,length(unique(legs_mne_crevise))),flag_mne(l+1:end)]; % cortical vs. subcortical flag

% Time Course Plots with Regional Legends
if ~isempty(time_norm)
sel_times = find(time_est < time_norm/1000);                                % x-indices for which to compute norm                                 
else
sel_times = 1:length(time_est);
end
[X_sim, n_sim, leg_sim, cflag_sim, X_greedy, n_greedy, leg_greedy, cflag_greedy, X_mne, n_mne, leg_mne, cflag_mne, ...
X_mne_all, n_mne_all, leg_mne_all, cflag_mne_all] = plot_timecourses_meg(Xvpl_plot, Xc_plot, Xls_plot_grcond, Xls_plot_mnecond_gr, ...
Xls_plot_mnecond_all, time_est, T, unique_rois, unique_cort_scort, legs_mne_condall, cort_scort_mne, col, thresh, fullcort, sel_times, sim_on); 

% Bar Graphs of Current Distributions
leg_sim{strcmp(leg_sim,'vpl')} = 'lt-vpl';                                  % fix legends to be intelligible
if fullcort == 0
    leg_greedy{strcmp(leg_greedy,'vpl')} = 'lt-vpl';                        % fix legends to be intelligible
    leg_mne{strcmp(leg_mne,'vpl')} = 'lt-vpl';                              % fix legends to be intelligible
    regions_list = [unique(rois(unique_patches)),strcat(regname,'-')];      %select cortical (based on sparse list) and all subcort regs
    cort_scort_flag = [ones(1,length(unique(rois(unique_patches)))) zeros(1, length(regname))];
else
    regions_list = [cort_rois, strcat(regname,'-')];                        %select cortical (based on shortlist from MNE) and all subcort regis
    cort_scort_flag = [ones(1,length(cort_rois)), zeros(1, length(regname))];
end
leg_greedy = strcat(leg_greedy,'-');                                        %so that comparisons are possible with updated regions list
currdistrn_rois(regions_list, n_sim, leg_sim, 'sim', n_greedy, leg_greedy, 'scsp',...
        cort_scort_flag, thresh2, fullcort, strcat(regname,'-'),1);         %bar graph for greedy estimates
currdistrn_rois(regions_list, n_sim, leg_sim, 'sim', n_mne, leg_mne, 'mne-sel', ...
        cort_scort_flag, thresh2, fullcort, strcat(regname,'-'),0);         %bar graph for mne estimates in selected regions
currdistrn_rois(regions_list, n_sim, leg_sim, 'sim', n_mne_all, leg_mne_all, 'mne', ...
        cort_scort_flag, thresh2, fullcort, strcat(regname,'-'),0);         %bar graph for all mne estimates

% Angles Across Hierarchical Stages
if fullcort == 0
plotangles_finalstage(lk, clkmax, Gtheta, c_curr_str, hipsurf_ind1, all_regions); % Check Angles for All Inverse Test Cases
end
diary off;