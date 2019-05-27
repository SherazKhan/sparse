%% Paths and Parameters
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

%Cortical Fwd Solutions
ico_num = 2;                                                                % cortical patch size
cfwd = mne_read_forward_solution(strcat(fwd_path, prefix, '-meg-all-fixed-fwd.fif'));% Source file for cortical forward solutions
load(strcat(fwd_path,prefix,'_cort_meg_',chtype,'_',gradmag,'_SVDs.mat'));  % Cortical forward solutions
load(strcat(fwd_path,'curr_den.mat'),'cort_cd');                            % Read in cortical current density
%Subcortical Fwd Solutions
srcspace_case = 2014;                                                       %'2013_lvpl_manual'; %2014_LT_submask1, 2015_LT_Submask5, 2016_LT_Submask2
sregname = 'lt';                                                            % Plot for left thalamus
load(strcat(fwd_path,'curr_den.mat'),'th_cd');                              % Read in VPL current density
load(strcat(fwd_path,prefix,'_ico3_forangles_meg_all_',gradmag,'.mat'),'Cw','meg_sel');% Read in whitener and MEG channels 
% whitener diagonal, preserves order + scort divisions to equal ico3 cort
p_sing = 6; nthresh = 95;                                                   % SVD computation parameters

%Fitting Parameters
nummodes = 1; nm = logical(nummodes == 1);                                  % Specify how many modes for the subcortical simulation
fittype = 'cort';                                                           % Specify fit type for all output files saved
SNR_fit = 9;                                                                % Used to determine lambda_sq = 1/SNR - previously set to 1

%Field Map Plotting Parameters
phys_meas_fname = strcat(subjdir, '/meas/',prefix,'_SomSens_SingleNoTask_offl_SSP.fif'); % For plotting field maps
virt_meas_fname = strcat(subjdir, '/meas/',prefix,'_virt_meg.fif');         % For plotting field maps
helmet_fname = '/autofs/eris/purdongp/users/pavitrak/sourceloc_scripts/ScSPIGH/fieldmaps_angles/Elekta_VirtMag304.mat';
proj_approach = 'cart_sph';                                                 % Choose 'cart_sph' or 'mdimsc
num_contours = 0;                                                           % Choose # of contours to display - was 50

%% Simulate Activity in a Subcortical Region
% Obtain Eigendecomposition of Prewhitened Fwd Solution for SrcSpace Cases Specified Above
if srcspace_case == 2014
  load(strcat(fwd_path, '2014-15_Vectorview_System/whiten_',gradmag, '/SM04_subc_meg_all_1g25m_SVDs.mat'),'sdiv_fwd');
  sregsdiv = 1;                                                             % Submask 1
end

% Identify the Division Index Corresponding to VPL Thalamus
for i = 1:length(sdiv_fwd) 
    all_regions{i} = sdiv_fwd(i).reg_name;                                  % All Region Names
end
ind = find(strcmp(sregname, all_regions));                                  % Identify all region names
sdiv = ind(sregsdiv); clear ind;                                            % Identify subdivision index amongst sdiv_fwd divisions

% Simulate Field Map in Most Significant Eigenmode
sp = 1*nm + sdiv_fwd(sdiv).numwhite_modes*(1-nm);                           % Relevant Number of Modes
Us = sdiv_fwd(sdiv).whiteUs(:,1:sp);                                        % Fwd Mode(s) of Interest
Ss = sdiv_fwd(sdiv).whiteS(1:sp,1:sp);                                      % Fwd Spectral Norm(s) of Interest
Gs = sdiv_fwd(sdiv).currstr*Us*Ss;                                          % Prewhitened Subcortical Fwd: Reduced Dim + Scaled by Currstr [Am]
if sp == 1
    xs = 1;                                                                 % Source Currents: Mode Space - 1 unit of currstr of this division (No Units)
else
    xs = 1./sqrt(diag(Ss)*sum(diag(Ss)));                                   % Source Currents: Mode Space - xs propto SingVals (No Units)
end
ys = Gs*xs;                                                                 % Data for Field Map : Prewhitened (No Units)
N = size(ys,1);                                                             % # Sensors

%% Resultant Amplitude of Simulated Activity in Subcortical Region
Vs = sdiv_fwd(sdiv).whiteVs(:,1:sp);                                        % Right eigenvector - num_dipoles X 1
Xs_plot.mode = xs*sdiv_fwd(sdiv).currstr;                                   % Mode Amplitudes [nAm]
Xs_plot.dip_allorient =  Vs*Xs_plot.mode;                                   % Dipole Amplitudes [Am]
Xs_plot.dip_or(:,:,1) = Xs_plot.dip_allorient(1:3:end,:);                   % Orientation 1
Xs_plot.dip_or(:,:,2) = Xs_plot.dip_allorient(2:3:end,:);                   % Orientation 2
Xs_plot.dip_or(:,:,3) = Xs_plot.dip_allorient(3:3:end,:);                   % Orientation 3
[Xs_plot.dip_res, Xs_plot.dip_resor] = vec_sum(Xs_plot.dip_or);             % Net simulated resultant activity for region [Am]
Xs_source = Xs_plot.dip_res;                                                % Resultant Source Current [Am]

%% Collate Forward Solutions of All Cortical Patches or Full Brain
cortp = p(ico_num);                                                         % ICO characteristics
cpatch = ico(ico_num).patch; L = length(cpatch); Gcfull = [];               % Initialization 
c_currstr = cort_cd*mean([cpatch(:).A]);                                    % mean cortical current strength across patches
for j = 1:L
    Gc = c_currstr*cpatch(j).U(:,1:cortp)*cpatch(j).S(1:cortp,1:cortp);     % Prewhitened Cortical Fwd: Reduced Dim + Scaled by Currstr [Am]
    patch_exc(j) = cpatch(j).exclude;                                       % Exclusion flag for medial wall patches
    if ~patch_exc(j)
        Gcfull = [Gcfull Gc];                                               % Use only the G's that are not on medial wall
    end
end
Gfit = Gcfull; lkmax = size(Gfit,2);                                        % Number of cortical modes to consider for fit
 
%% Fit Activity to All Cortical Patches
xfit = source_estimates('MNE', eye(N), Gfit, N, ys, 1, size(ys,2), SNR_fit, [], ones(size(Gfit,2)));
xfit_scaled = xfit*c_currstr;
Xc_plotmri = srcspace_activity_mri(prefix, L, xfit_scaled, cortp, cfwd, ico_num, source, ico, mri_path, patch_exc, lkmax, fittype, Xs_source);
scaled_patchres = [Xc_plotmri(:).dip_res_sc];                               % Mode to Dip, account orient, compute patch-wise resultants, scale by Xs_source
Xc_plotsurf = srcspace_activity_surf(prefix, L, xfit_scaled, cortp, cfwd, ico_num, source, ico, patch_exc, lkmax, fittype, Xs_source);
figure, set(gcf,'color','white'); plot(scaled_patchres); 
disp('run the following in terminal for current distribution plots')
disp('freeview -v $SUBJECTS_DIR/$SUBJECT/mri/brain.mgz SM04_ico_3_cortfit_mask.mgz $SUBJECTS_DIR/$SUBJECT/mri/svolume_decomps/140620_Thesis_SrcSpace/lt_submask1.mgz');

%% For Field Maps: DeWhiten to Physical Units
Gs_no_whiten = Cw^(-1)*Gs;                                                  % Dewhitening Gs implicitly dewhitens ys used for field maps
Gfit_no_whiten = Cw^(-1)*Gfit(:,1:lkmax);                                   % Dewhitening Gfit implicitly dewhitens ycfit used for field maps

%% Store For Field Maps (3D Physical Helmet)
% Field Map settings
data_3D = fiff_read_evoked_all(phys_meas_fname);                            % Read data file to modify
[chnames, chsel] = channel_selector('meg', data_3D, chtype);                % Select channels relevant to meastype
data_3D = fiff_pick_channels_evoked(data_3D, chnames);                      % Read data specifically for these channels
% Subcortical Activity Map
s_fmap_fname = strcat('real_s_mode_',fittype,'3D.fif');                     % Field map file name
fieldmap3D(data_3D, Gs_no_whiten, s_fmap_fname, chsel, xs);                 % Write field map for mne_analyze view      
% Cortical Fit Map
fc_fmap_fname = strcat('real_fc_mode_',fittype,'3D.fif');                   % Field map file name
fieldmap3D(data_3D, Gfit_no_whiten, fc_fmap_fname, chsel, xfit(1:lkmax,:)); % Write field map for mne_analyze view      

%% Convert Field Maps from Physical Helmet to Virtual Magnetometer Array 
news_fmap_fname = strcat('virt_s_mode_',fittype,'3D.fif');                  % Field map file name
newfc_fmap_fname = strcat('virt_fc_mode_',fittype,'3D.fif');                % Field map file name
disp('run the following in terminal for forward solutions');
unix_cmd{1} = ['mne_map_data --from ', s_fmap_fname, ' --frommri $SUBJECTS_DIR/$SUBJECT/fwd/$SUBJECT"-trans.fif" --to ',...
                virt_meas_fname, ' --tomri $SUBJECTS_DIR/$SUBJECT/fwd/$SUBJECT"-trans.fif" --out ', news_fmap_fname];
unix_cmd{2} = ['mne_map_data --from ', fc_fmap_fname, ' --frommri $SUBJECTS_DIR/$SUBJECT/fwd/$SUBJECT"-trans.fif" --to ', ...
                virt_meas_fname, ' --tomri $SUBJECTS_DIR/$SUBJECT/fwd/$SUBJECT"-trans.fif" --out ', newfc_fmap_fname];
disp(char(unix_cmd(~cellfun(@isempty,unix_cmd))));
disp('after it completes: type exit at terminal');                          % Should be Fast
pause; clc;

%% Plot Field Maps (2D Virtual Magnetometer)
% Subcortical Activity Map
s_data_2D = fiff_read_evoked_all(news_fmap_fname);                          % Subcortical data in virtual mag system
s_data_2Dtoplot = s_data_2D.evoked.epochs(:,1);                             % Only 1 time point required
s_scale = max(abs(min(s_data_2Dtoplot)), abs(max(s_data_2Dtoplot)));  	    % Normalized subcortical field map colorscale
s_data_2Dtoplot_sc = s_data_2Dtoplot/s_scale;				    % Normalized subcortical field map
s_ax = fieldmap2D(helmet_fname, proj_approach, virt_meas_fname, s_data_2Dtoplot_sc, num_contours, []);  
set(gcf,'renderer','painters'); print(gcf,'-depsc','fmap_subcortical_source');% -depsc for color
% Cortical Fit Map
fc_data_2D = fiff_read_evoked_all(newfc_fmap_fname);                        % Cortical fit in virtual mag system
fc_data_2Dtoplot = fc_data_2D.evoked.epochs(:,1);                           % Only 1 time point required
c_scale = max(abs(min(fc_data_2Dtoplot)), abs(max(fc_data_2Dtoplot)));	    % Normalized cortical field map colorscale
fc_data_2Dtoplot_sc = fc_data_2Dtoplot/c_scale;				    % Normalized cortical field map	
fc_ax = fieldmap2D(helmet_fname, proj_approach, virt_meas_fname, fc_data_2Dtoplot_sc, num_contours, s_ax);  
set(gcf,'renderer','painters'); print(gcf,'-depsc','fmap_cortical_fit');    % -depsc for color

%% Compute Metrics to Report
numer = norm(s_data_2Dtoplot - fc_data_2Dtoplot).^2;                        % numerator for GOF  - MSE
denom = norm(s_data_2Dtoplot).^2;                                           % denominator for GOF - MSE
gof = 1 - sqrt(numer/denom);                                                % GOF b/w Subcortical Source & Best Fit Cortical Field Maps
pang = max(subspacea(Gs, Gfit))*180/pi;                                     % principal angle for this specific combination
disp(sdiv_fwd(sdiv).volume); disp(gof); disp(pang);			    % metrics to report