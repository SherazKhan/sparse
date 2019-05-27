clear; clc; close all;
addpath '/autofs/eris/purdongp/users/pavitrak/sourceloc_data/AEP/cat_004/meas/abr_plr_runall'; 
addpath '/autofs/eris/purdongp/users/pavitrak/sourceloc_data/AEP/cat_004/meas/abr_plr_calibration'; 
covs = 1;

%% Average Measurement Files
raw_fname = 'eyesopen_1_raw.fif';
raw = fiff_setup_read_raw(raw_fname);
abrmeas_fname = 'abr_runall-HFABR-bc-ave.fif'; 
%abrmeas_fname = 'abr_runall-ABR-bc-ave.fif'; 
abrmeas = fiff_read_evoked_all(abrmeas_fname);
abrY = abrmeas.evoked.epochs;
mlrmeas_fname = 'abr_runall-MLR-bc-ave.fif';
mlrmeas = fiff_read_evoked_all(mlrmeas_fname);
mlrY = mlrmeas.evoked.epochs;
abrLeff = double(abrmeas.evoked.nave); 
mlrLeff = double(mlrmeas.evoked.nave); 

%% Pick MEG Channels
want_meg = true;
want_eeg = false;
want_stim = false;
include = {};
exclude = raw.info.bads;
megpicks = fiff_pick_types(raw.info,want_meg,want_eeg,want_stim,include,exclude);
megch = {raw.info.ch_names{megpicks}};

%% Pick EEG Channels
want_meg = false;
want_eeg = true;
want_stim = false;
include = {};
exclude = raw.info.bads;
eegpicks = fiff_pick_types(raw.info,want_meg,want_eeg,want_stim,include,exclude);
eegch = {raw.info.ch_names{eegpicks}};

%% Choose Channels for Given Region
selpath = '/autofs/homes/008/pavitrak/.mne/';
selname = strcat('cat004.txt');
A = importdata(selname);
%str = {'Right-occipital','Left-occipital'};
%str = {'Right-temporal','Left-temporal'};
%str = {'Right-temporal','Left-temporal','Right-parietal','Left-parietal'};
str = {'Right-temporal','Left-temporal','Right-occipital','Left-occipital','Right-frontal','Left-frontal'};
cnt = 1;
for ss = 1:length(str)
    selstr = str{ss}; 
    row = strmatch(selstr,A);
    importdata(selname)
    selch = A{row}(length(selstr)+2:end);
    ind = strfind(selch,'S');
    numch = (ind-1)/9;
    for i = cnt:cnt+numch-1
        ii = i-cnt+1;
        st = ii + 8*(ii-1);
        chname{i} = selch(st:st+7);
        chname{i} = strcat(chname{i}(1:3),chname{i}(5:end));
    end
    cnt = length(chname) + 1;
end

%% Read Relevant Parts of Measurement Files for MEG/EEG
if covs == 0
    selmegch = megch(ismember(megch,chname)); 
else
    selmegch = megch; 
end
abrmeas_meg = fiff_pick_channels_evoked(abrmeas,selmegch,exclude);
abrmeas_eeg = fiff_pick_channels_evoked(abrmeas,eegch,exclude);
mlrmeas_meg = fiff_pick_channels_evoked(mlrmeas,selmegch,exclude);
mlrmeas_eeg = fiff_pick_channels_evoked(mlrmeas,eegch,exclude);
if covs == 0
    coil_type = [abrmeas_meg.info.chs(:).coil_type];%selmegch
else
    coil_type = [abrmeas_meg.info.chs(:).coil_type];
end
cij = ismember(coil_type,3024);

%% MEG/EEG Time Courses
Tma = abrmeas_meg.evoked.times;
Yma = abrmeas_meg.evoked.epochs;
Tea = abrmeas_eeg.evoked.times;
Yea = abrmeas_eeg.evoked.epochs;
Tmm = mlrmeas_meg.evoked.times;
Ymm = mlrmeas_meg.evoked.epochs;
Tem = mlrmeas_eeg.evoked.times;
Yem = mlrmeas_eeg.evoked.epochs;

%% Plot Time Courses
figure(1), set(gcf,'color','white');  
subplot(3,1,1), plot(Tea*1000, mean(abs(Yea))*10^(6)'); 
title('Mean Rectified ABR Measurement','FontSize',14);
xlabel('Time (Milliseconds)','FontSize',14); ylabel('EEG Amplitude (uV)','FontSize',14);
subplot(3,1,2), plot(Tma*1000, mean(abs(Yma(cij,:)))*10^(15)'); 
title('MEG Magnetometer','FontSize',14);
xlabel('Time (Milliseconds)','FontSize',14); ylabel('Amplitude (fT)','FontSize',14);
subplot(3,1,3), plot(Tma*1000, mean(abs(Yma(~cij,:)))*10^(13)'); 
title('MEG Gradiometer','FontSize',14);
xlabel('Time (Milliseconds)','FontSize',14); ylabel('Amplitude (fT/cm)','FontSize',14);

figure(2), set(gcf,'color','white'); 
subplot(3,1,1), plot(Tem*1000, Yem*10^(6)'); 
title('MLR Measurement','FontSize',14); 
xlabel('Time (Milliseconds)','FontSize',14); ylabel('EEG Amplitude (uV)','FontSize',14);
subplot(3,1,2), plot(Tmm*1000, Ymm(cij,:)*10^(15)'); 
title('MEG Magnetometer','FontSize',14);
xlabel('Time (Milliseconds)','FontSize',14); ylabel('Amplitude (fT)','FontSize',14)
subplot(3,1,3), plot(Tmm*1000, Ymm(~cij,:)*10^(13)');
title('MEG Gradiometer','FontSize',14)
xlabel('Time (Milliseconds)','Fontsize',14); ylabel('Amplitude (fT/cm)','FontSize',14)

%% Covariances Filenames 
mlr = str2double(input('MLR (1) or ABR (0)')); 
megcov = str2double(input('MEG (1) or EEG (0)'));
mlrbccov_fname= 'eyesopen_MLR-bc-cov.fif';  %eyesopen_12_MLR12-bc-cov.fif
abrbccov_fname= 'eyesopen_HFABR-bc-cov.fif';%eyesopen_12_HFABR12-bc-cov.fif
%abrbccov_fname= 'eyesopen_ABR-bc-cov.fif'; %eyesopen_12_ABR12-bc-cov.fif
rawssp_fname  = 'emptyroom_raw.fif';        
mlrssp_fname  = 'emptyroom_MLR_SSP.fif';    
abrssp_fname  = 'emptyroom_HFABR_SSP.fif';
if mlr == 1
    cov_fname = mlrbccov_fname;
    ssp_fname = mlrssp_fname;
    Leff = mlrLeff;
else
    cov_fname = abrbccov_fname;
    ssp_fname = abrssp_fname;
    Leff = abrLeff;
end

%% Read in Noise Covariances
noise = mne_read_noise_cov(cov_fname);
if megcov == 1
    noise = mne_pick_channels_cov(noise,megch);
else
    noise = mne_pick_channels_cov(noise,eegch);
end
C0 = noise.data; Cn = C0/Leff;

%% Project MEG SSP onto MEG Covs
if megcov == 1
    [fid,dir] = fiff_open(mlrssp_fname); 
    projdata = fiff_read_proj(fid,dir); fclose(fid);
    for ii = 1:length(projdata)
        projdata(ii).active = 1; 
    end
    P = mne_make_projector(projdata,megch);
else
    if mlr
    eegprojdata = mlrmeas_eeg.info.projs(end).data.data;
    else
    eegprojdata = abrmeas_eeg.info.projs(end).data.data;
    end
    a = eegpicks -306; a= [a(1:54) a(55:64)-4];
    eegprojdata = eegprojdata(a); 
    P = eye(length(eegch)) - eegprojdata'*eegprojdata;
end    
Cp = P*Cn*P;

%% Load Covariances
if megcov == 1
diagC = diag(Cp);
Ekm = 0.20;                   
Ekg = 0.05;
Cl = Cp + Ekm*diag(mean(diagC(cij))*double(cij)) + Ekg*diag(mean(diagC(~cij))*double(~cij));
else
diagC = diag(Cp);
Eke = 0.05;
Cl = Cp + Eke*diag(mean(diagC)*ones(length(eegch),1));
end
[V,D] = eig(Cl);
lambda2p = diag(D);
for ii = 1:length(lambda2p)
	if isreal(lambda2p(ii)) && lambda2p(ii) > 0
        %can try to use threshold instead of 0, 1.5*10^(-31)
		lambda2p_f(ii) = lambda2p(ii);
	end
end
lambda2p = lambda2p_f;

%% Plots of Covariances - without and with SSP, and Eigenvalues
if megcov == 1
 figure(5), set(gcf, 'color','white'); imagesc(Cn(cij, cij));%sqrt(abs(
 title('Magnetometer Covariance Structure - without SSP','FontSize',14); 
 colorbar; set(gca,'FontSize',14); cmag = caxis;
 xlabel('Magnetometer Index','FontSize',14); ylabel('Magnetometer Index','FontSize',14);
 
 figure(6), set(gcf, 'color','white'); imagesc(Cn(~cij, ~cij));
 title('Gradiometer Covariance Structure - without SSP','FontSize',14); 
 colorbar; set(gca,'FontSize',14); cgrad = caxis;
 xlabel('Gradiometer Index','FontSize',14); ylabel('Gradiometer Index','FontSize',14); 
 
 figure(7), set(gcf, 'color','white'); imagesc(Cp(cij, cij));
 title('Magnetometer Covariance Structure - with SSP','FontSize',14); 
 caxis(cmag); colorbar; set(gca,'FontSize',14); 
 xlabel('Magnetometer Index','FontSize',14); ylabel('Magnetometer Index','FontSize',14); 
 
 figure(8), set(gcf, 'color','white'); imagesc(Cp(~cij, ~cij));
 title('Gradiometer Covariance Structure - with SSP','FontSize',14);  colorbar;
 caxis(cgrad); colorbar; set(gca,'FontSize',14);
 xlabel('Gradiometer Index','FontSize',14); ylabel('Gradiometer Index','FontSize',14); 
 
 figure(9), set(gcf,'color','white'); s = sort(lambda2p); semilogy(s);
 title('Eigenvalues of MEG Covariance Matrix','FontSize',14)
 set(gca,'FontSize',14);
 ylabel('Eigenvalues','FontSize',14);  xlabel('Eigenvalue Index','Fontsize',14)
else
 figure(5), set(gcf,'color','white'); imagesc(Cn); 
 title('EEG Covariance Structure - without SSP','FontSize',14);
 ceeg = caxis; colorbar; set(gca,'FontSize',14);
 xlabel('EEG Channel Index','FontSize',14); ylabel('EEG Channel Index','FontSize',14); 
 
 figure(6), set(gcf,'color','white'); imagesc(Cp); 
 title('EEG Covariance Structure - with SSP','FontSize',14);
 caxis(ceeg); colorbar; set(gca,'FontSize',14);
 xlabel('EEG Channel Index','FontSize',14); ylabel('EEG Channel Index','FontSize',14);
 
 figure(7), set(gcf,'color','white'); s = sort(lambda2p); semilogy(s);
 title('Eigenvalues of EEG Covariance Matrix','FontSize',14)
 set(gca,'FontSize',14);
 ylabel('Eigenvalues','FontSize',14);  xlabel('Eigenvalue Index','FontSize',14)
end