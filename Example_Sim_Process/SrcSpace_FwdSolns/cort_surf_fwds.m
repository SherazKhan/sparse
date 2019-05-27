%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Load  Cortical Fwd Solns and Compute Associated SVDs
%Written by Pavitra Krishnaswamy 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Load Cortical Fwd Solns 
megfwd_fname = strcat(datapath,'/fwd/',prefix,'-meg-all-fixed-fwd.fif');    %Vectorview Forward solution
megfwd_mindist_fname = strcat(datapath,'/fwd/',prefix,'-meg-all-mindist-fixed-fwd.fif'); %FWd after excluding superficial patches
llabs = strcat(prefix,'-meg-all-mindistout-lh-label.txt');                  %labels of left patches to exclude
rlabs = strcat(prefix,'-meg-all-mindistout-rh-label.txt');                  %labels of right patches to exclude
fwd_meg = mne_read_forward_solution(megfwd_fname);                          %Hippocampal ico MEG fwd solution
fwd_meg_mindist = mne_read_forward_solution(megfwd_mindist_fname);          %Hippocampal ico mindist fwd solution
mindistout_labels(prefix, fwd_meg, fwd_meg_mindist,llabs, rlabs);           %Labels < mindist
crawG = fwd_meg.sol.data;                                                   %Raw G matrix to be split
meas_fname = strcat(datapath, '/meas/',prefix,'_vpl_somato_meg_sim',append);%measurement file

%% Load Whitener                                                        
load(meas_fname,'Cw');                                                      %Need Whitened SVD to calculate effective current strength!

%% Whiten Fwds
meg_sel = 1:1:306;                                                          %Selected MEG Channels
cwhiteG = Cw*crawG(meg_sel,:);                                              %Whitened cortical forward solution

%% Compute SVDs
csvd_name = strcat(datapath,'/fwd/',prefix,'_cort-',meastype,'-all_white_SVDs',append,'.mat'); 
p_sing = 6;                                                                 %#Singular values to compute per patch 
pt_s = 10;                                                                  %Perc. superficial dipoles tolerable
nthresh = 95;                                                               %NMRA threshold for #modes per patch
subdiv = 3;                                                                 %Run patch decomposition all through to ico-4
if exist(csvd_name,'file') ~= 2 
  %note fwd_meg is only used for its src structures - real measurement info is in crawG
  cpatch_decomp(prefix, fwd_meg, cwhiteG, nthresh, p_sing, pt_s, subdiv, csvd_name, meastype);
end
load(csvd_name,'ico','p'); eigmodes = p; clear p;                           %whitened SVD file
