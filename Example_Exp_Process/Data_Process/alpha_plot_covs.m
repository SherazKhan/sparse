clear; clc; close all;
addpath '/autofs/eris/purdongp/users/pavitrak/sourceloc_data/AEP/cat_004/meas/abr_plr_runall';
addpath '/autofs/eris/purdongp/users/pavitrak/sourceloc_data/AEP/cat_004/meas/abr_plr_calibration';

%% Read in MEG Noise Covariances
raw_fname = 'emptyroom_raw.fif';
raw = fiff_setup_read_raw(raw_fname);
ssp_fname = 'emptyroom_EYEC_SSP.fif';
cov_fname = 'emptyroom_EYEC-bc-cov.fif';
want_meg = true;
want_eeg = false;
want_stim = false;
include = {};
exclude = raw.info.bads;
megpicks = fiff_pick_types(raw.info,want_meg,want_eeg,want_stim,include,exclude);
megch = {raw.info.ch_names{megpicks}};
noise = mne_read_noise_cov(cov_fname);
noise = mne_pick_channels_cov(noise,megch);
Leff = 1; C0 = noise.data; Cn = C0/Leff;

%% Project SSP onto Covs
[fid,dir] = fiff_open(ssp_fname); 
projdata = fiff_read_proj(fid,dir); fclose(fid);
for ii = 1:length(projdata)
    projdata(ii).active = 1; 
end
P = mne_make_projector(projdata,megch);   
Cp = P*Cn*P;

%% Load Covariances
diagC = diag(Cp);
Ekm = 0.20;                   
Ekg = 0.05;
coil_type = [raw.info.chs(megpicks).coil_type];
cij = ismember(coil_type,3024);
Cl = Cp + Ekm*diag(mean(diagC(cij))*double(cij)) + Ekg*diag(mean(diagC(~cij))*double(~cij));
[V,D] = eig(Cl);
lambda2p = diag(D);
for i = 1:length(lambda2p)
	if isreal(lambda2p(i)) == 0 || lambda2p(i) < 0
		lambda2p(i) = 0;
    end
end

%% Plots of Covariances - without and with SSP, and Eigenvalues
figure(1), set(gcf, 'color','white'); imagesc(Cn(cij, cij));%sqrt(abs(
title('Magnetometer Covariance Structure - without SSP','FontSize',14); 
colorbar; set(gca,'FontSize',14);
xlabel('Magnetometer Index','FontSize',14); ylabel('Magnetometer Index','FontSize',14); 
figure(2), set(gcf, 'color','white'); imagesc(Cn(~cij, ~cij));
title('Gradiometer Covariance Structure - without SSP','FontSize',14);  
colorbar; set(gca,'FontSize',14);
xlabel('Gradiometer Index','FontSize',14); ylabel('Gradiometer Index','FontSize',14); 
figure(3), set(gcf, 'color','white'); imagesc(Cp(cij, cij));
title('Magnetometer Covariance Structure - with SSP','FontSize',14); 
colorbar; set(gca,'FontSize',14);
xlabel('Magnetometer Index','FontSize',14); ylabel('Magnetometer Index','FontSize',14); 
figure(4), set(gcf, 'color','white'); imagesc(Cp(~cij, ~cij));
title('Gradiometer Covariance Structure - with SSP','FontSize',14);  
colorbar; set(gca,'FontSize',14);
xlabel('Gradiometer Index','FontSize',14); ylabel('Gradiometer Index','FontSize',14); 
figure(5), set(gcf,'color','white'); s = sort(lambda2p); semilogy(s);
title('Eigenvalues of MEG Covariance Matrix','FontSize',14); 
ylabel('Eigenvalues','FontSize',14);  xlabel('Eigenvalue Index','FontSize',14)