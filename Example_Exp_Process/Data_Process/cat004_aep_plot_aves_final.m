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
%str = {'Right-temporal','Left-temporal','Right-occipital','Left-occipital'};
str = {'Vertex','Left-temporal','Right-temporal','Left-parietal','Right-parietal',...
        'Left-occipital','Right-occipital','Left-frontal','Right-frontal'};
cnt = 1;
for ss = 1:length(str)
    selstr = str{ss}; 
    row = strmatch(selstr,A);
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
xlabel('Time (Milliseconds)','Fontsize',14); ylabel('Amplitude (fT/cm)','FontSize',14);

%% Plot EEG Raw Time Course
figure, set(gcf,'color','white'); 
plot(Tea*1000, Yea*10^6');
ylabel('EEG amplitude (uV)'); xlabel('Time (Milliseconds)');