%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Load  Cortical Fwd Solns and Compute Associated SVDs
%Written by Pavitra Krishnaswamy 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Load Whitener, Channel Selection                                                      
meas_fname = strcat(datapath,'/meas/',prefix, '-',meastype, '-', datatype, append, '_measinfo.mat'); %whitener, covariance, channels sel, proj
load(meas_fname,'meg_sel','megCw','megP','eeg_sel','eegCw','eegP','exclude');%Need Whitened SVD to calculate effective current strength!
if strcmp(meastype,'meg') 
	sel = meg_sel;  Cw = megCw;     P = megP; 
elseif strcmp(meastype,'eeg')
    sel = eeg_sel;  Cw = eegCw;     P = eegP; 
end

%% Load Cortical Fwd Solns 
meg_fwd_fname = strcat(prefix, '-meg_runavg-all-fixed-fwd.fif');
eeg_fwd_fname = strcat(prefix, '-eeg-all-fixed-fwd.fif');
meg_fwd_mindist_fname = strcat(prefix, '-meg_runavg-all-mindist-fixed-fwd.fif');
eeg_fwd_mindist_fname = strcat(prefix, '-eeg-all-mindist-fixed-fwd.fif');
fwd_meg = mne_read_forward_solution(meg_fwd_fname, {}, {}, {}, exclude);
fwd_eeg = mne_read_forward_solution(eeg_fwd_fname, {}, {}, {}, exclude);
fwd_meg_mindist = mne_read_forward_solution(meg_fwd_mindist_fname, {}, {}, {}, exclude);
fwd_eeg_mindist = mne_read_forward_solution(eeg_fwd_mindist_fname, {}, {}, {}, exclude);
llabs = strcat(prefix,'-meg-all-mindistout-lh-label.txt');                  %labels of left patches to exclude
rlabs = strcat(prefix,'-meg-all-mindistout-rh-label.txt');                  %labels of right patches to exclude
mindistout_labels(prefix, fwd_meg, fwd_meg_mindist,llabs, rlabs);           %Labels < mindist

%% Whiten Fwds
cmeg_whiteG = megCw*(megP*fwd_meg.sol.data);                                % Apply SSP projection and whiten meg fwd solution
ceeg_whiteG = eegCw*(eegP*fwd_eeg.sol.data);                                % Apply Average Reference Projection and whiten eeg fwd solution
cwhiteG = cat(1, cmeg_whiteG,ceeg_whiteG);                                  % concatenate meg and eeg fwds

%% Compute SVDs
csvd_name = strcat(datapath,'/fwd/',prefix,'_cort-',meastype,'_all-white_SVDs',append,'.mat');
p_sing = 6;                                                                 %#Singular values to compute per patch 
pt_s = 10;                                                                  %Perc. superficial dipoles tolerable
nthresh = 95;                                                               %NMRA threshold for #modes per patch
subdiv = 3;                                                                 %Run patch decomposition all through to ico-4
if exist(csvd_name,'file') ~= 2 
  %note fwd_meg is only used for its src structures - real measurement info is in crawG
  cpatch_decomp(prefix, fwd_meg, cwhiteG, nthresh, p_sing, pt_s, subdiv, csvd_name, meastype);
end
load(csvd_name,'ico','p'); eigmodes = p; clear p;                           %whitened SVD file