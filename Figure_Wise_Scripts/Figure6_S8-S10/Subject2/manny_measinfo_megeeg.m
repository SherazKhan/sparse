function manny_measinfo_megeeg(prefix, datatype, meastype, datapath, plot_on, append, saveon)

%% Specify Measurement, Covariance and Projector File Names
meas_fname = strcat(datapath,'/meas/abr_plr_runall/', 'abr_runall-',append(2:end),'-bc-ave.fif');
cov_fname = strcat(datapath, '/meas/abr_plr_calibration/', 'eyesopen',append,'-bc-cov.fif');
ssp_fname = strcat(datapath, '/meas/abr_plr_calibration/', 'emptyroom',append,'_SSP.fif');
if saveon
diary(strcat(prefix, '_', datatype, '_', meastype, append,'_matlog'));
savename = strcat(datapath,'/meas/',prefix, '-',meastype, '-', datatype, append, '_measinfo.mat');
end

%%  Make Channel Selector (Exclude Bad Channels)
% Select Good MEG Channels
meas = fiff_read_evoked_all(meas_fname);                                    % meas gives bad chs
want_meg = true;  want_eeg = false;  want_stim = false;                     % include only meg channels
include = {};  ear_lobes = {'EEG097' 'EEG098'};  exclude = [meas.info.bads ear_lobes]; % include and exclude channels
meg_sel = fiff_pick_types(meas.info, want_meg, want_eeg, want_stim, include, exclude); % selected meg channels
meg_chs = meas.info.ch_names(meg_sel);                                      % names of megchs included

% Select Good EEG Channels
want_meg = false;  want_eeg = true;  want_stim = false;                     % include only eeg channels
eeg_sel = fiff_pick_types(meas.info, want_meg, want_eeg, want_stim, include, exclude);
eeg_chs = meas.info.ch_names(eeg_sel);                                      % names of eegchs included
eeg_sel = [eeg_sel(1:58) eeg_sel(59:length(eeg_chs))-4]-(min(eeg_sel)-(max(meg_sel)+1)); % update EEG Sel for later Use
%eeg_sel = [eeg_sel(1:58) eeg_sel(59:length(eeg_chs))-4];                   % update EEG Sel for later Use - not same as above for this subject

% Joint Measurement
meas = fiff_pick_channels_evoked(meas, [meg_chs eeg_chs]);                  % meas for good channels only

%%  Read and Pre-Process Measurement
Y = meas.evoked.epochs(:,:);                                                % ERP measurement
Fs = meas.info.sfreq;                                                       % sampling frequency
t0 = 101;                                                                   % start of file to start of stim
delay = 0;                                                                  % start of stim to start of erp
Y = Y(:,t0+delay:end);                                                      % selected ERP window
%ceil(Fs*(t0+delay)):floor(Fs*(t0+delay+tw))                                % tw = length(ERP response in sec)
N = size(Y,1);                                                              % number channels
t = size(Y,2);                                                              % number time points
t1 = 1; t2 = t;

%%  Create STC Structure Fields
stc(1).tmin = delay;                                                        % seconds
stc(2).tmin = stc(1).tmin;                                                  % seconds
stc(1).tstep = 1/Fs;                                                        % seconds
stc(2).tstep = stc(1).tstep;                                                % seconds
time = linspace(0,t-1,t)*stc(1).tstep + delay;                              % time vector
%time = linspace(ceil(Fs*stc(1).tmin),ceil(Fs*stc(1).tmin)+(t-1),t)*stc(1).tstep;

%%  Eyes Open Noise Covariance
noise = mne_read_noise_cov(cov_fname);                                      % read in noise covariance
Leff = double(meas.evoked.nave);                                            % number of effective epochs

% MEG
meg_noise = mne_pick_channels_cov(noise,meg_chs,eeg_chs);                   % read noise cov for selected ch
megC0 = meg_noise.data;                                                     % noise covariance
megCn = megC0/Leff;                                                         % scale by # epochs

% EEG
eeg_noise = mne_pick_channels_cov(noise,eeg_chs,meg_chs);                   % read noise cov for selected ch
eegC0 = eeg_noise.data;                                                     % noise covariance
eegCn = eegC0/Leff;                                                         % scale by # epochs

%%  Read Empty Room SSP, Make SSP Projector and Apply to MEG Noise Cov.
[fid,dir] = fiff_open(ssp_fname);                                           % projector filename
meg_projdata = fiff_read_proj(fid, dir);                                    % read projector matrix
fclose(fid);
for ii = 1:length(meg_projdata)
    meg_projdata(ii).active = 1;                                            % activate projectors
end
megP = mne_make_projector(meg_projdata, meg_chs);                           % make projector
megCp = megP*megCn*megP;                                                    % apply projector

%%  Make EEG Average Reference Projector and Apply to EEG Noise Cov.
N_eeg = length(eeg_chs);
eegP = eye(N_eeg) - ones(N_eeg)/N_eeg;                                      % projector
eegCp = eegP*eegCn*eegP;                                                    % projection
eegCp = diag(diag(eegCp));                                                  % diagonalize EEG noise covariance

%%  Load or Regularize MEG Noise Covariance
N_meg = length(meg_chs);
coil_type = [meas.info.chs(:).coil_type];                                   % coils differentiate mags vs grads
cij = ismember(coil_type, [3022 3024]);                                     % mag indexing by cij and grad indexing by ~cij
%cij = ismember(coil_type,3024);                                            %mag indexing by cij
cij = cij(1:N_meg);                                                         % so that cij only takes MEG channels into account
Ekm = 0.50;                                                                 % magnetometer loading
Ekg = 0.05;                                                                 % gradiometer loading
diag_megCp = diag(megCp);
megCl = megCp + Ekm*diag(mean(diag_megCp(cij))*double(cij)) + Ekg*diag(mean(diag_megCp(~cij))*double(~cij));

%%  Load or Regularize EEG Noise Covariance
Eke = 0.05;                                                                 % EEG loading
eegCl = eegCp + Eke*diag(mean(diag(eegCp))*ones(N_eeg,1));

%%  Derive EigenDecomposition of MEG Noise Covariance and Compute MEG Whitener
V = zeros(N_meg);                                                           % initialize eigenvectors
lambda2_meg = zeros(1, N_meg);                                              % initialize eigenvalues
lambda2_meg_inv = zeros(N_meg);                                             % initialize inv of eigenvalues

% Magnetometer block
[V_mag, D_mag] = eig(megCl(cij,cij));                                       % eigendecomposition for magnetometers
V(cij,cij) = V_mag;                                                         % eigenvectors for magnetometers
lambda2_meg(cij) = diag(D_mag);                                             % eigenvalues for magnetometers
m_i = strfind({meg_projdata.desc}, 'axial');                                % details of projnmag setting
mag_n_eigs = sum([m_i{:}]);                                                 % projnmag setting
mag_min_eigs = [];
if ~isempty(mag_n_eigs)
    mag_ord_eigs = sort(diag(D_mag), 'descend');                            % sort eigenvalues for magnetometers
    mag_min_eigs = mag_ord_eigs(end:-1:end-(mag_n_eigs-1));                 % lowest projnmag eigenvalues for magnetometers
end

% Gradiometer block
[V_grad, D_grad] = eig(megCl(~cij,~cij));                                   % eigendecomposition for gradiometers
V(~cij,~cij) = V_grad;                                                      % eigenvectors for gradiometers
lambda2_meg(~cij) = diag(D_grad);                                           % eigenvalues for gradiometers
g_i = strfind({meg_projdata.desc}, 'planar');                               % details of projngrad setting
grad_n_eigs = sum([g_i{:}]);                                                % projngrad setting
grad_min_eigs = [];
if ~isempty(grad_n_eigs)
    grad_ord_eigs = sort(diag(D_grad), 'descend');                          % sort eigenvalues for gradiometers
    grad_min_eigs = grad_ord_eigs(end:-1:end-(grad_n_eigs-1));              % lowest projngrad eigenvalues for gradiometers
end

% Form whitening operator (megCw)
for ii = 1:length(lambda2_meg)
    if isreal(lambda2_meg(ii)) && lambda2_meg(ii) > 0 && ~ismember(lambda2_meg(ii), [mag_min_eigs(:)' grad_min_eigs(:)'])
        lambda2_meg_inv(ii,ii) = 1./sqrt(lambda2_meg(ii));
    end
end
megCw = lambda2_meg_inv*V';                                                 % meg whitener

%%  Derive EigenDecomposition of EEG Noise Covariance and Compute EEG Whitener
[V_eeg, D_eeg] = eig(eegCl);                                                % eigendecomposition of eeg noise cov
lambda2_eeg = diag(D_eeg);                                                  % eigenvalues
lambda2_eeg_inv = zeros(N_eeg);                                             % initialize inv of eigenvalues
lamthreshe = 0;  % lamthreshe = min(lambda2_eeg);                           % lowest real eigenvalues for EEG is 0
for ii = 1:length(lambda2_eeg)
    if isreal(lambda2_eeg(ii)) && lambda2_eeg(ii) > 0 && ~ismember(lambda2_eeg(ii), lamthreshe)
        lambda2_eeg_inv(ii,ii) = 1./sqrt(lambda2_eeg(ii));
    end
end
eegCw = lambda2_eeg_inv*V_eeg';                                             % eeg whitener

%%  Plot Covs and Eigenvalues
if plot_on
    figure, set(gcf,'color','white'), imagesc(megCn(~cij,~cij)); colorbar;
    figure, set(gcf,'color','white'), imagesc(megCn(cij,cij)); colorbar;
    figure, set(gcf,'color','white'), imagesc(megCp(~cij,~cij)); colorbar;
    figure, set(gcf,'color','white'), imagesc(megCp(cij,cij));  colorbar;
    figure, set(gcf,'color','white'), imagesc(megCl(~cij,~cij)); colorbar;
    figure, set(gcf,'color','white'), imagesc(megCl(cij,cij));  colorbar;
    figure, set(gcf,'color','white'), imagesc(eegCn); colorbar;
    figure, set(gcf,'color','white'), imagesc(eegCp); colorbar;
    figure, set(gcf,'color','white'), imagesc(eegCl); colorbar;
    if strcmp(append, '_HFABR')
        caxis([0 0.015]*10^-14);
    end
    figure, set(gcf,'color','white'), semilogy(sort(real(lambda2_meg(~ismember(lambda2_meg, [mag_min_eigs(:)' grad_min_eigs(:)'])))));
    figure, set(gcf,'color','white'), semilogy(sort(real(lambda2_eeg(~ismember(lambda2_eeg, lamthreshe)))));
    figure, set(gcf,'color','white'); imagesc(megCw(cij,cij)); colorbar;
    figure, set(gcf,'color','white'); imagesc(megCw(~cij,~cij)); colorbar;
    figure, set(gcf,'color','white'); imagesc(eegCw); colorbar;
    pause
    close all
end

%% Whiten Data and Forward Solution
megYraw = Y(1:N_meg,:);                                                     % raw meg measurement 
megYp = megP*megYraw;                                                       % Apply SSP projection to meg measurement
eegYraw = Y(N_meg+1:end,:);                                                 % raw eeg measurement
eegYp = eegP*eegYraw;                                                       % Apply Average Reference projection to eeg measurement
megYw = megCw*megYp;                                                        % whiten meg measurement
eegYw = eegCw*eegYp;                                                        % whiten eeg measurement
Yraw = cat(1,megYraw,eegYraw);                                              % concatenate raw meg and eeg
Y = cat(1, megYw, eegYw);                                                   % concatenate meg and eeg measurements
C = eye(N_meg+N_eeg);                                                       % noise cov after whitening
scalarX = 56;                                                               % for computing resultants during visualization

%% Save Important Variables and Call Inverse Solutions
%if strcmp(mnetype, 'stable')                                
%	meg_sel = 1:1:size(megCw,1);                                            %need to exclude only earlobes from fwd, other bad chs already accounted for
%end
if saveon
save(savename, 'megYraw','megYp','eegYraw','eegYp','megYw','eegYw','Yraw',...
    'exclude','cij', 'meg_sel', 'eeg_sel', 'Ekg', 'Ekm', 'Eke', ...
    'megC0', 'megCn', 'megP', 'megCl', 'megCw', 'eegC0', 'eegCn', 'eegP', 'eegCl', 'eegCw',...
    'Y', 'N','Fs', 't1', 't2', 't', 'Leff', 'stc', 'time','C','scalarX','-v7.3');
diary off;
end
           
end