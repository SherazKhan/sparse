function svd_for_fwds(prefix, meastype, datatype, ico_num, append, varargin)

%% Set Paths and Load Forward Solutions
% General Paths and Parameters
clc; close all;
add_allpaths;
meas_fname = strcat(datapath, '/meas/', prefix, '-',meastype, '-', datatype, append, '_measinfo.mat'); %meas info
if nargin == 6
    whiten_fname = varargin{1};                                             %whitener file name
else
    whiten_fname = meas_fname;                                              %whitener file name
end

% Cortical Info
meg_fwd_fname = strcat(prefix, '-meg_runavg-all-fixed-fwd.fif');            %cortical forward solution
eeg_fwd_fname = strcat(prefix, '-eeg-all-fixed-fwd.fif');
meg_fwd_mindist_fname = strcat(prefix, '-meg_runavg-all-mindist-fixed-fwd.fif');%cortical forward solution with 5mm sources excluded
eeg_fwd_mindist_fname = strcat(prefix, '-eeg-all-mindist-fixed-fwd.fif');
csvds_fname = strcat(datapath, '/fwd/', prefix,'_cort_',meastype,'_plr-prewhite','_SVDs', append,'.mat');
subdiv_stop = 4;  redo_csvd = 0;                                            %run cortical ico decomposition till this subdiv

% Subcortical Info
add_on_name = []; %add_on_name = '_bs-full';
ssvds_fname = strcat(datapath, '/fwd/', prefix,'_subc_',meastype,'_plr-prewhite','_SVDs',add_on_name, append,'.mat');
volumes_fname = strcat(datapath, '/fwd/',prefix,'_subc_',meastype,'_plr-prewhite_volumes',add_on_name, append,'.mat');
redo_ssvd = 1; 

% Angles Files
angs_fname = strcat(datapath,'/fwd/',prefix,'_forangles_',meastype,'_plr-prewhite',add_on_name, append,'.mat');
redo_angs = 1;                                                              %whether to redo SVD even if the files exist

%% SVD Computation for Fwd Matrix on Selected Channels
% Parameters
p_sing = 6;                                                                 % #singular values to compute per patch
pt_s = 100;                                                                 % perc. superficial dipoles tolerable
nthresh = 95;                                                               % NMRA threshold for #modes per patch
exclude_regs = {'lh','rh'};                                                 %do not include hippocampal volumes
thresh = 0.90;                                                              %Require > 90% of Patch Intersect with Chosen Regions
load(whiten_fname,'megCw','megP','eegCw','eegP','exclude','meg_sel','eeg_sel');%Whitener (prewhitened SVD to calculate effective current strength), other info

% Cortical SVDs
if exist(csvds_fname,'file') ~= 2 || redo_csvd == 1
    fwd_meg = mne_read_forward_solution(meg_fwd_fname, {}, {}, {}, exclude);%read fwd soln   
    fwd_eeg = mne_read_forward_solution(eeg_fwd_fname, {}, {}, {}, exclude);
    fwd_meg_mindist = mne_read_forward_solution(meg_fwd_mindist_fname, {}, {}, {}, exclude); %read fwd mindist soln
    fwd_eeg_mindist = mne_read_forward_solution(eeg_fwd_mindist_fname, {}, {}, {}, exclude);
    mindistout_labels(prefix, fwd_meg, fwd_meg_mindist);                    %determine superficial labels
    meg_whiteG = megCw*(megP*fwd_meg.sol.data);                             % Apply SSP projection and whiten meg fwd solution
    eeg_whiteG = eegCw*(eegP*fwd_eeg.sol.data);                             % Apply Average Reference Projection and whiten eeg fwd solution
    whiteG = cat(1, meg_whiteG, eeg_whiteG);                                % concatenate meg and eeg fwds
    cpatch_decomp(prefix, fwd_meg, whiteG, nthresh, p_sing, pt_s, subdiv_stop, csvds_fname, meastype);
    for i = 1:subdiv_stop
        excl_medial_labels(datapath, prefix, ico_num, csvds_fname, thresh); %Exclude Medial Patches
    end
end
load(csvds_fname,'source','ico','p');                                       %eigendecompositions of cortical patch fwds
cpatch = ico(ico_num).patch; clear ico;                                     %reduced order cortical patch fwds
lpatch_ind = length(source(ico_num).hemisph(1).pinfo);                      %indices marking left patches  

% Subcortical SVDs
if strcmp(append,'_HFABR')
  if exist(ssvds_fname,'file') ~= 2 || redo_ssvd == 1 
    [sregname, sdiv_fwd, best_raw_sdiv_modes, best_white_sdiv_modes]...     %eigendecompositions of subcortical volume fwds
    = svolume_decomp_sdiv(prefix, meastype, datapath, volumes_fname, megCw*megP, eegCw*eegP, meg_sel, eeg_sel, nthresh, p_sing);
    for i = 1:length(sdiv_fwd)  
      chk_exclude(i) = sum(strcmp(sdiv_fwd(i).reg_name, exclude_regs));     %exclude regions not of interest
    end
    sdiv_fwd = sdiv_fwd(chk_exclude == 0);                                  %select regions which are not flagged for exclusion
    for j = 1:length(exclude_regs)
        sregname = sregname(~strcmp(sregname, exclude_regs{j}));            %select region names which are not flagged for exclusion
    end
    save(ssvds_fname,'sregname','sdiv_fwd','best_raw_sdiv_modes','best_white_sdiv_modes','-v7.3');
  end
  load(ssvds_fname,'sdiv_fwd');
end

%% Reduce Order of Eigendecomposition (NMRA Criterion) to Match Our Inverse Analysis
if strcmp(append,'_HFABR')
  % Cortical Eigendecompositions
  cortp = p(ico_num); clear p; 
  for i = 1:length(cpatch)
    cpatch(i).name = i;                                                     %patch index
    if i < lpatch_ind
        cpatch(i).hemi = 'l';                                               %left hemisphere
    else
        cpatch(i).hemi = 'r';                                               %right hemisphere
    end
    cpatch(i).cp = cortp;                                                   %number of modes for NMRA
    if cpatch(i).exclude == 0
      cpatch(i).U = cpatch(i).U(:, 1:cpatch(i).cp);                         %left eigenvector
      cpatch(i).S = cpatch(i).S(1:cpatch(i).cp, 1:cpatch(i).cp);            %eigenvalues
      cpatch(i).V = cpatch(i).V(:, 1:cpatch(i).cp);                         %right eigenvector
    else
      cpatch(i).U = []; cpatch(i).S = []; cpatch(i).V = [];                 %dont store anything as medial patches
    end
  end

  % Subcortical Eigendecompositions
  for i = 1:length(sdiv_fwd)
    sdiv(i).name = sdiv_fwd(i).reg_name;                                    %name of region
    if strfind(sdiv_fwd(i).reg_name(1),'l') == 1
        sdiv(i).hemi = 'l';                                                 %left hemisphere
    elseif strfind(sdiv_fwd(i).reg_name(1),'r') == 1
        sdiv(i).hemi = 'r';                                                 %right hemisphere
    end
    sdiv(i).sp = sdiv_fwd(i).numwhite_modes;                                %number of modes for NMRA
    sdiv(i).U = sdiv_fwd(i).whiteUs(:, 1:sdiv(i).sp);                       %left eigenvector
    sdiv(i).S = sdiv_fwd(i).whiteS(1:sdiv(i).sp, 1:sdiv(i).sp);             %eigenvalues
    sdiv(i).V = sdiv_fwd(i).whiteVs(:, 1:sdiv(i).sp);                       %right eigenvector
  end

  if exist(angs_fname,'file') ~= 2 || redo_angs == 1
  save(angs_fname,'cpatch','lpatch_ind','sdiv','p_sing', 'nthresh', 'pt_s','megCw','eegCw',...
     'meg_sel','eeg_sel','csvds_fname', 'ico_num', 'ssvds_fname');
  end
end

end