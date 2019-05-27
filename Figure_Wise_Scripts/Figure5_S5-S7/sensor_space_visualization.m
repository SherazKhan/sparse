%% Paths and Parameters
% Subject and Simulation File Parameters
clear; clc; close all;
addpath(genpath('/autofs/eris/purdongp/users/pavitrak/sourceloc_scripts/ScSPIGH'));
prefix = 'SM04'; meastype = 'meg'; datatype = 'SEP';                        % subject name
subjdir = strcat('/autofs/eris/purdongp/users/pavitrak/sourceloc_data/',datatype,'/',prefix); % subject data folder
sim_fname = strcat(subjdir,'/meas/',prefix,'_vpl_somato_meg_sim_overlap.mat'); % simulation file
load(sim_fname,'time','Y');                                                 % load simulation data

% Field Map Generation and Plotting Parameters
virt_meas_fname = strcat(subjdir, '/meas/',prefix,'_virt_meg.fif');         % For plotting field maps
helmet_fname = '/autofs/eris/purdongp/users/pavitrak/sourceloc_scripts/ScSPIGH/fieldmaps_angles/Elekta_VirtMag304.mat';
proj_approach = 'cart_sph';                                                 % Choose 'cart_sph' or 'mdimsc
num_contours = 0;                                                           % Choose # of contours to display - was 50
simfield_fname = strcat(subjdir,'/meas/',prefix,'_SomSens_Simulation.fif'); % Simulation Evoked FIF File
seltime = 84;                                                          	    % snapshot sensor fields at 84 msecs
	
%% Read Data (Evoked Fiff) File to Modify
chtype = 'all';                                                             % Consider both magnetometers and gradiometers
phys_meas_fname = strcat(subjdir, '/meas/',prefix,'_SomSens_SingleNoTask_offl_SSP.fif'); % For plotting field maps
phys_data = fiff_read_evoked_all(phys_meas_fname);                          % Read data file to modify
[chnames, chsel] = channel_selector('meg', phys_data, chtype);              % Select channels relevant to meastype
phys_data = fiff_pick_channels_evoked(phys_data, chnames);                  % Read data specifically for these channels

%% Generate a Simulated Data File and Write Modified Evoked File
simdata.info = phys_data.info;                                              % copy information and channel fields
simdata.info.sfreq = 1/(time(2)-time(1));                                   % update sampling frequency
simdata.evoked = phys_data.evoked;                                          % copy evoked fields
%time = time(1:3:end); Y = Y(:,1:3:end);                                    % resample time and mesurement matrices
simdata.evoked.first = 1; %time(1)*1000;                                        % first recording time
simdata.evoked.last = 792; %time(end)*1000;                                       % last recording time
simdata.evoked.times = time*1000;                                           % ERP recording times 
simdata.evoked.epochs = Y(chsel,:);                                         % only include y for selected channels
fiff_write_evoked(simfield_fname, simdata);                                 % write simulated field map file

%% Convert Field Maps from Physical Helmet to Virtual Magnetometer Array 
simvirt_fmap_fname = strcat(pwd, '/', 'sim_fieldmap_virt3D.fif');           % Field map file name
disp('run the following in terminal for forward solutions');
unix_cmd{1} = ['mne_map_data --from ', simfield_fname, ' --frommri $SUBJECTS_DIR/$SUBJECT/fwd/$SUBJECT"-trans.fif" --to ',...
                virt_meas_fname, ' --tomri $SUBJECTS_DIR/$SUBJECT/fwd/$SUBJECT"-trans.fif" --out ', simvirt_fmap_fname];
disp(char(unix_cmd(~cellfun(@isempty,unix_cmd))));
disp('after it completes: type exit at terminal'); pause; clc;

%% Plot Field Maps (2D Virtual Magnetometer)
sim_data_2D = fiff_read_evoked_all(simvirt_fmap_fname);                     % Subcortical data in virtual mag system
[~, closest_ind] = min(abs(simdata.evoked.times - seltime));                % nearest index to selected time point
sim_data_2Dtoplot = sim_data_2D.evoked.epochs(:,closest_ind);               % Selected time point
sim_ax = fieldmap2D(helmet_fname, proj_approach, virt_meas_fname, sim_data_2Dtoplot, num_contours, []);  
set(gcf,'renderer','painters'); print(gcf,'-depsc',strcat('simulated_fmap'));% -depsc for color