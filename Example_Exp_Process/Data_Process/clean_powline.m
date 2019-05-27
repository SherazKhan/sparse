%% Set Paths and Parameters 
clear; clc; close all;
addpath(genpath('/autofs/eris/purdongp/users/pavitrak/matlab/')); startup(1);
datapath = '/autofs/eris/purdongp/users/pavitrak/sourceloc_data/AEP/cat_004/meas/preproc_rawdata/';
addpath(datapath);
evoked = 1; raw1 = 0; ploton = 0; 
evoked_fnames = {'ABR_run1','ABR_run2','ABR_run3','ABR_run4','ABR_run5'};
non_evoked_fnames = {'eyesopen_1','eyesopen_2','eyesclosed_1','eyesclosed_2','emptyroom','emptyroom_stimon'};

for jj = 1:5 
  if evoked == 1
    if raw1 == 0
        meas_fname0 = strcat(evoked_fnames{jj},'_raw.fif');
    else
        meas_fname0 = strcat(evoked_fnames{jj},'_raw-1.fif');
    end
  else
    meas_fname0 = strcat(non_evoked_fnames{jj},'_raw.fif');
  end
  meas_fname = strcat(datapath, meas_fname0);
  cax = [-320 -280]; NW = 3; K = 2*NW-1; movingwin = [4 4]; specres = NW*2/4

  %% Read in All Data Channels
  raw = fiff_setup_read_raw(meas_fname);
  Fs = raw.info.sfreq; 
  params = struct('tapers',[NW K], 'pad', -1, 'Fs', Fs, 'fpass', [0 Fs/2],'err',0);
  Y = fiff_read_raw_segment(raw);
  N = size(Y,1); t = size(Y,2); 
  chan = 1:1:380; %don't apply to MISC and STI channels
  newdata = zeros(size(Y));

  for cc = 1:length(chan)
  %% Read Data
  t1 = 1; t2 = t; ch = chan(cc); 
  disp(strcat('run  ', num2str(jj),' ...channel  ', num2str(ch)));
  data_fs(:,1) =  (t1:1:t2)/Fs;            
  data_fs(:,2) =  detrend(Y(ch,t1:t2));
  Estim.Raw = data_fs(:,2);
  
  %% Plot Raw Spectrogram
  [Spect.Raw.PSD,Spect.Raw.Times,Spect.Raw.Freqs] = mtspecgram(Estim.Raw', params, movingwin, ploton);
  if ploton == 1
      set(gcf,'color','white'); title('Raw Data'); cax = caxis; colorbar; 
  end
 
  %% Form Notch Filters and Apply Serially to Data
  y = Estim.Raw; 
  for i = 1:30
    wo = 60*i/(Fs/2);  bw = 1/(Fs/2); %wo/35;
    [bnotch,anotch] = iirnotch(wo,bw);
    a = anotch; b = bnotch; 
    %fvtool(b,a); pause; close; 
    y0 = filter(b,a,y);
    y = y0;
  end
  Estim.Resid = y; newdata(ch,:) = Estim.Resid; 
  Estim.Harm = Estim.Raw - Estim.Resid; %size(Estim.Resid)
  
  %% Plot Harmonic Estimate and Cleaned Spectrogram
  [Spect.Harm.PSD,Spect.Harm.Times,Spect.Harm.Freqs] = mtspecgram(Estim.Harm', params, movingwin, ploton);
  if ploton == 1
      set(gcf,'color','white'); title('Harmonic Estimate'); caxis(cax); colorbar;  
  end
  [Spect.Resid.PSD,Spect.Resid.Times,Spect.Resid.Freqs] = mtspecgram(Estim.Resid', params, movingwin, ploton);
  if ploton == 1
      set(gcf,'color','white'); title('Power Line Noise Removed'); caxis(cax); colorbar; 
  end
  %spect_powline(ch) = Spect; estim_powline(ch) = Estim;
  
  if ploton == 1
      pause; close all; 
  end
  clear data_fs y y0 Estim spect_harm stimes_harm Spect
  end

  %% Write FIF File for Processing After Power Line Removal
  newdata(381:385,:) = Y(381:385,:);
  infile = meas_fname;
  outfile = strcat(pwd, '/powline/',meas_fname0);
  mne_modify_fiff(infile, outfile, newdata);
  %save('test.mat','spect_powline','estim_powline');
  clear meas_fname meas_fname0 raw params Y N t t1 t2 chan infile outfile newdata 
end