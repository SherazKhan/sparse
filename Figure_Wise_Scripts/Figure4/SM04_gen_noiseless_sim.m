%% Paths and Parameters
clear; clc; close all;
% General Parameters
addpath(genpath('/autofs/eris/purdongp/users/pavitrak/sourceloc_scripts/ScSPIGH'));
prefix = 'SM04'; meastype = 'meg'; datatype = 'Illustrations';
subjdir = strcat('/autofs/eris/purdongp/users/pavitrak/sourceloc_data/',datatype,'/',prefix);
g = 1; m = 25; gradmag = strcat(num2str(g),'g',num2str(m),'m'); chtype = 'all';
savepath = strcat(pwd,'/2016_Results_MatsFigs/');
% Simulation Parameters
load(strcat(savepath, prefix,'-',meastype,'-',datatype,'_sinv.mat'));
load(csvds_fname,'ico','p'); ico = ico(cdiv); eigmodes = p(cdiv);  clear p; % forward solutions: patch decompositions
ico_num = cdiv; lk = unique(patch_est); cpatch = ico.patch;                 % Read in Data
load(strcat(subjdir, '/fwd/curr_den.mat'), 'cort_cd');  sim_p = 1;          % Relevant Current Strength and # Modes for Simulation
load(ssvds_fname,'sdiv_fwd','sregname');                                    % subcortical SVDs
% Read Source Spaces and Forward Solutions
csrc = mne_read_source_spaces(strcat(subjdir, '/src/', prefix, '-ico-',num2str(cdiv), 'p-src.fif'));% cort source space
cfwd = mne_read_forward_solution(strcat(subjdir, '/fwd/', prefix, '-meg-all-fixed-fwd.fif'));% cort fwd
hsrc = mne_read_source_spaces(strcat(subjdir, '/src/', prefix, '_hipsurf-ico-1p-src.fif'));  % hip source space
hfwd = mne_read_forward_solution(strcat(subjdir, '/fwd/', prefix, '_hipsurf-meg-all-fixed-fwd.fif'));% hip fwd

%% Simulate Activity in a Cortical Region
c_curr_str = cort_cd*mean([cpatch(:).A]);                                   % average current strength for any cortical patch
for i = 1:length(lk)
Uc{i} = cpatch(lk(i)).U(:,1:sim_p);                                         % Fwd Mode of Interest
Sc{i} = cpatch(lk(i)).S(1:sim_p,1:sim_p);                                   % Fwd Spectral Norm of Interest
Gc{i} = Uc{i}*Sc{i};                                                        % Reduced Cortical Forward Solution (not scaled by currstr)
xc{i} = 1*c_curr_str;                                                       % Source Currents: Mode Space [Am]
Y(:,i) = Gc{i}*xc{i};                                                       % Data for Field Map [physical units]
end
%%%%Put in Xsource resultant calculation for cortical patches here%%%%%%%%
count = size(Y,2);

for i = count+1:count+length(sdiv_fwd)
%% Simulate Activity in a Subcortical Region
all_regions{i} = sdiv_fwd(i-count).reg_name;                                % All Region Names
sdiv = sdiv_fwd(i-count);                                                   % Sdiv Structure for this Region
Us{i} = sdiv.whiteUs(:,1:sim_p);                                            % Fwd Mode of Interest
Ss{i} = sdiv.whiteS(1:sim_p,1:sim_p);                                       % Fwd Spectral Norm of Interest
Gs{i} = Us{i}*Ss{i};                                                        % Reduced Subcortical Forward Solution (not scaled by currstr)
xs{i} = 1*sdiv.currstr;                                                     % Source Currents: Mode Space [Am]
Y(:,i) = Gs{i}*xs{i};                                                       % Data for Field Map [physical units]

%% Resultant Amplitude of Simulated Activity in Subcortical Region
Vs = sdiv.whiteVs(:,1:sim_p);                                               % Right eigenvector - num_dipoles X 1
Xs_plot(i).mode = xs{i};                                                    % Mode Amplitudesn [Am]
Xs_plot(i).dip_allorient =  Vs*Xs_plot(i).mode;                             % Dipole Amplitudes [Am]
if isempty(strfind(all_regions{i},'hipsurf'))                               % Not hippocampus surface
  Xs_plot(i).orient = [1 1 1];                                              % Orientation vector
  Xs_plot(i).dip_or(:,:,1) = Xs_plot(i).dip_allorient(1:3:end,:);           % Orientation 1
  Xs_plot(i).dip_or(:,:,2) = Xs_plot(i).dip_allorient(2:3:end,:);           % Orientation 2
  Xs_plot(i).dip_or(:,:,3) = Xs_plot(i).dip_allorient(3:3:end,:);           % Orientation 3
else                                                                        % Hippocampal surface
  Xs_plot(i).dip = Xs_plot(i).dip_allorient;                                % Dipole current estimate
  if isempty(strfind(all_regions{i},'rh'));                                 % left orientation vectors
  [~,~,hlk] = intersect([hsrc(1).pinfo{[sdiv.index]}],hfwd.src(1).vertno,'legacy'); 
  Xs_plot(i).orient = hsrc(1).nn(hlk,:);                   
  %rr{i} = hsrc(1).rr(hlk,:);
  else                                                                      % right orientation vectors
  [~,~,hrk] = intersect([hsrc(2).pinfo{[sdiv.index]}],hfwd.src(2).vertno,'legacy'); 
  Xs_plot(i).orient = hsrc(2).nn(hrk,:);                   
  %rr{i} = hsrc(2).rr(hrk,:);
  end
  nn(:,:,1) = Xs_plot(i).orient; %nn = repmat(nn,[1,1,T]);                  % Orientation vector
  Xs_plot(i).dip_or(:,:,1) = squeeze(nn(:,1,:)).*Xs_plot(i).dip;            % Orientation 1
  Xs_plot(i).dip_or(:,:,2) = squeeze(nn(:,2,:)).*Xs_plot(i).dip;            % Orientation 2
  Xs_plot(i).dip_or(:,:,3) = squeeze(nn(:,3,:)).*Xs_plot(i).dip;            % Orientation 3
  clear nn
end
[Xs_plot(i).dip_res, Xs_plot(i).dip_resor] = vec_sum(Xs_plot(i).dip_or);    % Net simulated resultant activity for region [Am]
Xs_source(i) = Xs_plot(i).dip_res;                                          % Resultant Source Current [Am]
end

%% Store Noiseless Simulated Measurements
N = size(Y,1);                                                              % # Sensors
t1 = 1; t2 = 1; t = 1;                                                      % # Time points
C = eye(N);                                                                 % Noise covariance for MNE - as all simulations in whitened modes
sim_fname = strcat(savepath,prefix,'-', meastype, '-',datatype,'_noiseless_sim.mat');
save(sim_fname,'sim_p','Uc','Sc','Gc','xc','all_regions','Us','Ss','Gs','xs','Xs_plot','Xs_source',...
    'Y','C','N','t1','t2','t');
save(strcat(savepath,prefix,'-',meastype,'-',datatype,'_sinv.mat'),'sim_p','Uc','Sc','Gc','xc',...
    'all_regions','Us','Ss','Gs','xs','Y','Xs_plot','Xs_source','N','t1','t2','t','-append');