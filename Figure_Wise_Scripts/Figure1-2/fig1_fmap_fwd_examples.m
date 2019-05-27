%% Paths and Parameters
% General Parameters
clear; clc; close all;
addpath(genpath('/autofs/eris/purdongp/users/pavitrak/sourceloc_scripts/ScSPIGH'));
prefix = 'SM04'; meastype = 'meg'; datatype = 'Illustrations';
subjdir = strcat('/autofs/eris/purdongp/users/pavitrak/sourceloc_data/',datatype, '/',prefix);
chtype = 'all';                                                             % Consider both magnetometers and gradiometers
fwd_path = strcat(subjdir, '/fwd/');                                        % Forward solutions are all stored here

% Forward Solution
ico_num = 2;                                                                % size of cortical patches - visual choice
cfwd = mne_read_forward_solution(strcat(fwd_path, prefix, '-meg-all-fixed-fwd.fif'));% source file for cortical forward solutions
g = 1; m = 25;                                                              % Scaling ratios for gradiom vs. magnetom
gradmag = strcat(num2str(g),'g',num2str(m),'m');                            % Naming string
csvds_fname = strcat(fwd_path,prefix,'_cort_meg_',chtype,'_',gradmag,'_SVDs.mat');
thresh = 0.90;                                                              %Require > 90% of Patch Intersect with Chosen Regions
%excl_medial_labels(subjdir, prefix, ico_num, csvds_fname, thresh);          %Exclude Medial Patches
load(csvds_fname);                                                          % Cortical forward solutions
add_on_name = []; %add_on_name = '_bs-full';                                % Reduced or full brainstem
load(strcat(fwd_path,prefix,'_subc_ico3_meg_',chtype,'_',gradmag,'_SVDs',add_on_name,'.mat'));  %subcortical forward solutions - divided to equal ico3 cort
load(strcat(fwd_path,prefix,'_ico3_forangles_meg_all_',gradmag,'.mat'),'Cw','meg_sel'); 	% Read in whitener (diagonal, preserves order) and MEG channels
load(strcat(fwd_path,'curr_den.mat'),'cort_cd');                            % Cortical Current Density
csrc =  mne_read_source_spaces(strcat(subjdir, '/src/',prefix, '-ico-',num2str(ico_num),'p-src.fif'));

% Field Map Settings
helmet_fname = '/autofs/eris/purdongp/users/pavitrak/sourceloc_scripts/ScSPIGH/fieldmaps_angles/Elekta_VirtMag304.mat';
proj_approach = 'cart_sph';                                                 % Choose 'cart_sph' or 'mdimsc
num_contours = 0;                                                           % Choose # of contours to display
phys_meas_fname = strcat(subjdir, '/meas/',prefix,'_SomSens_SingleNoTask_offl_SSP.fif'); % For plotting field maps
virt_meas_fname = strcat(subjdir, '/meas/',prefix,'_virt_meg.fif');         % For plotting field maps

%% Simulate Activity in a Subcortical Region
if isempty(add_on_name)
    sregname = 'bsred'; sregsdiv = 9;%17;                                   % Relevant Region and Sdiv - Reduced Brainstem
else
    sregname = 'bs';    sregsdiv = 85;                                      % Relevant Region and Sdiv - Full Brainstem 
end
for i = 1:length(sdiv_fwd )
    all_regions{i} = sdiv_fwd (i).reg_name;                                 % All Region Names
end
ind = find(strcmp(sregname, all_regions)); sdiv = ind(sregsdiv); clear ind; % Identify subdivision index;
sp = 1;                                                                     % Pick top mode only
Us = sdiv_fwd(sdiv).whiteUs(:,1:sp);                                        % Fwd Mode(s) of Interest
Ss = sdiv_fwd(sdiv).whiteS(1:sp,1:sp);                                      % Fwd Spectral Norm(s) of Interest
Gs = sdiv_fwd(sdiv).currstr*Us*Ss;                                          % Prewhitened Subcortical Fwd: Reduced Dim + Scaled by Currstr [Am]
xs = 1;                                                                     % Source Currents: Mode Space - 1 unit of currstr of this division (No Units)
ys = Gs*xs;                                                                 % Data for Field Map : Prewhitened (No Units)

%% Simulate Activity in a Cortical Patch
cpatch = ico(ico_num).patch; L = length(cpatch); Gcfull = [];               % Initialization 
c_currstr = cort_cd*mean([cpatch(:).A]);                                    % mean cortical current strength across patches
patch_exc = 1;
while patch_exc
    patchno = 100;%str2num(input('pick a cortical patch j:  '));            % pick a patch j
    patch_exc = cpatch(patchno).exclude;                                    % Exclusion flag for medial wall patches
end
clabel_no = patch_to_labloc(csrc, cfwd, patchno)
cp = 1; %cortp = p(ico_num);                                                % Pick top mode only
Gc = c_currstr*cpatch(patchno).U(:,1)*cpatch(patchno).S(1,1);               % Prewhitened Cortical Fwd: Reduced Dim + Scaled by Currstr [Am]
xc = 1;                                                                     % Source Currents: Mode Space - 1 unit of currstr of this division (No Units)
yc = Gc*xc;                                                                 % Data for Field Map : Prewhitened (No Units)

%% For Field Maps: DeWhiten to Physical Units
Gs_no_whiten = Cw^(-1)*Gs;                                                  % Dewhitening Gs implicitly dewhitens ys used for field maps
Gc_no_whiten = Cw^(-1)*Gc;                                                  % Dewhitening Gc implicitly dewhitens yc used for field maps

%% Store For Field Maps (3D Physical Helmet)
% Field Map settings
data_3D = fiff_read_evoked_all(phys_meas_fname);                            % Read data file to modify
[chnames, chsel] = channel_selector('meg', data_3D, chtype);                % Select channels relevant to meastype
data_3D = fiff_pick_channels_evoked(data_3D, chnames);                      % Read data specifically for these channels
% Subcortical Activity Map
s_fmap_fname = strcat('real_s_mode_3D.fif');                                % Field map file name
fieldmap3D(data_3D, Gs_no_whiten, s_fmap_fname, chsel, xs);                 % Write field map for mne_analyze view      
% Cortical Activity Map
c_fmap_fname = strcat('real_c_mode_3D.fif');                                % Field map file name
fieldmap3D(data_3D, Gc_no_whiten, c_fmap_fname, chsel, xc);                 % Write field map for mne_analyze view      

%% Convert Field Maps from Physical Helmet to Virtual Magnetometer Array 
news_fmap_fname = strcat('virt_s_mode_3D.fif');                             % Field map file name
newc_fmap_fname = strcat('virt_c_mode_3D.fif');                             % Field map file name
disp('run the following in terminal for forward solutions');
unix_cmd{1} = ['mne_map_data --from ', s_fmap_fname, ' --frommri $SUBJECTS_DIR/$SUBJECT/fwd/$SUBJECT"-trans.fif" --to ',...
                virt_meas_fname, ' --tomri $SUBJECTS_DIR/$SUBJECT/fwd/$SUBJECT"-trans.fif" --out ', news_fmap_fname];
unix_cmd{2} = ['mne_map_data --from ', c_fmap_fname, ' --frommri $SUBJECTS_DIR/$SUBJECT/fwd/$SUBJECT"-trans.fif" --to ', ...
                virt_meas_fname, ' --tomri $SUBJECTS_DIR/$SUBJECT/fwd/$SUBJECT"-trans.fif" --out ', newc_fmap_fname];
disp(char(unix_cmd(~cellfun(@isempty,unix_cmd))));
disp('after it completes: type exit at terminal'); pause; clc;

%% Plot Field Maps (2D Virtual Magnetometer)
% Subcortical Activity Map
s_data_2D = fiff_read_evoked_all(news_fmap_fname);                          % Subcortical data in virtual mag system
s_data_2Dtoplot = s_data_2D.evoked.epochs(:,1);                             % Only 1 time point required
s_ax = fieldmap2D(helmet_fname, proj_approach, virt_meas_fname, s_data_2Dtoplot, num_contours, []);  
set(gcf,'renderer','painters'); print(gcf,'-depsc',strcat('fig1_fmap_subcortical_example',add_on_name));% -depsc for color

% Cortical Fit Map
c_data_2D = fiff_read_evoked_all(newc_fmap_fname);                          % Cortical fit in virtual mag system
c_data_2Dtoplot = c_data_2D.evoked.epochs(:,1);                             % Only 1 time point required
c_ax = fieldmap2D(helmet_fname, proj_approach, virt_meas_fname, c_data_2Dtoplot, num_contours, []);                         
set(gcf,'renderer','painters'); print(gcf,'-depsc','fig1_fmap_cortical_example');% -depsc for color