%% Paths and Filenames
clear; clc; close all; 
addpath '/autofs/eris/purdongp/users/pavitrak/sourceloc_data/AEP/cat_004/meas/abr_powline/abr_calibration';
ldfname = 'left_ear_tube_delay.fif';    
rdfname = 'right_ear_tube_delay.fif';   
lst = 1;        len = 100000;   lamp = 0.4; %lst = 20;  len = 60
rst = 120000;   ren = 200000;   ramp = 0.4; %rst = 0;   ren = 40;
lr = 1;
if lr == 1
    dfname = ldfname; st = lst; en = len; amp = lamp;
else
    dfname = rdfname; st = rst; en = ren; amp = ramp;  
end

%% Read 2 MISC Channels from Delay File
raw = fiff_setup_read_raw(strcat(dfname));
%meas = fiff_read_raw_segment_times(raw, st, en);
%will have to crosscheck that en < (raw.last_samp - raw.first_samp)/raw.info.sfreq
meas = fiff_read_raw_segment(raw, raw.first_samp, raw.last_samp);
ind_src = find(strncmp('MISC006',raw.info.ch_names, 7)==1);
ind_ear = find(strncmp('MISC004',raw.info.ch_names, 7)==1);
snd_src = meas(ind_src,:);
snd_ear = meas(ind_ear,:);

%% Visualize Delay
figure, set(gcf,'color','white');
snd_src = detrend(snd_src(st:en)); 
snd_ear = detrend(snd_ear(st:en));
subplot(2,1,1), plot(snd_src); ylabel('SRC');
subplot(2,1,2), plot(snd_ear); ylabel('EAR');
linkaxes;

%% Compute Delay (Milliseconds)
a  = snd_src - snd_ear;
z = zeros(size(a));
z(find(a>amp)) = 1;
figure, set(gcf,'color','white');
plot(a,'LineWidth',3);
hold on; plot(z,'r-');

figure, set(gcf,'color','white');
ind = find(z==1); 
dind = diff(ind); hist(dind,length(dind))
del = mean(dind(dind > 40 & dind < 60))/raw.info.sfreq*1000

save('delay.mat','del');