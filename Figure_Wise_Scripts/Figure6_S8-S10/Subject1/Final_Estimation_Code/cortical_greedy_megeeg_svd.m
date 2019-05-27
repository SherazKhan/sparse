%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SPIGH SOLUTION TO THE NEUROMAGNETIC INVERSE PROBLEM
%Localize Cortical Sources with Subspace Pursuit (Reduced State-Space Model)
%Written by Gabriel Obregon (obregon@nmr.mgh.harvard.edu)
%Modified and Adapted by Pavitra Krishnaswamy (pavitrak@nmr.mgh.harvard.edu)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Set Paths, Load Cortical Forward Solutions, Meas Files
% Basic Parameters
clear; clc; close all; 
prefix = 'cat_004'; meastype = 'meg_eeg'; datatype = 'AEP'; append = '_MLR';
diary(strcat(prefix, '_', meastype, '_cort_matlog',append));
%add_allpaths;
datapath = strcat('/autofs/eris/purdongp/users/pavitrak/sourceloc_data/',datatype,'/', prefix,'/');
addpath(genpath('/autofs/eris/purdongp/users/pavitrak/sourceloc_scripts/ScSPIGH/')); 
addpath('/autofs/eris/purdongp/users/pavitrak/sourceloc_scripts/figure_making/Figure7-8_S7-S9/');

% Specify Files with Meas, Noise and Fwd Info
measinfo_fname = strcat(datapath, 'meas/', prefix,'-',meastype,'-',datatype,append,'_measinfo.mat'); %meas, covariance, whitener
meg_fwd_fname = strcat(datapath, 'fwd/', prefix, '-meg_runavg-all-fixed-fwd.fif');            %MEG fwd solution file
eeg_fwd_fname = strcat(datapath, 'fwd/', prefix, '-eeg-all-fixed-fwd.fif');                   %EEG fwd solution file
csvds_fname = strcat(datapath, 'fwd/', prefix,'_cort_',meastype,'_plr-prewhite','_SVDs',append,'.mat'); %prewhitened SVDs
curr_den_fname = strcat(datapath, 'fwd/curr_den.mat');                      %for cortical current density

% Load Files with Meas, Noise and Fwd Info
load(measinfo_fname,'exclude','Y','C','N','t1','t2','stc','time','scalarX');%relevant measurement information
fwd_meg = mne_read_forward_solution(meg_fwd_fname, {}, {}, {}, exclude);    %MEG fwd solution
fwd_eeg = mne_read_forward_solution(eeg_fwd_fname, {}, {}, {}, exclude);    %EEG fwd solution
load(csvds_fname); eigmodes = p;  clear p;                                  %forward solutions: patch decompositions

%% Set Greedy Parameters
pt_s = 100;                                                                 %perc. superficial dipoles tolerable
alg = 'MNE'; %'dSPM'; %'sMAP-EM'                                            %if sMAP-EM do gradiometers only
if strcmp(alg,'sMAP-EM');
    stop = 7;
else
    stop = [];
end
subdiv_end = 3;                                                             %SPIGH to ICO-3 (could also run MNE to ICO-4 after)
SNR = 1;                                                                    %typical ERP case 3^2
L = 4; dub = 0;                                                             %target sparsity level and doubling heuristic
patch_est = cell(1,min(subdiv_end,4)); time_est = patch_est; 
Xdip = struct('patch',{},'alg',{});
ploton = 1; scale_currstr = 0;                                              %scale patch-wise if scale_currstr set to 1

%% Compute average worst-case mutual coherence b/w any randomly-selected patch and its nearest neighbors (u)
disp('mutual coherence with neighbors only');
u = cort_mutcoh_neighbors(csvds_fname,eigmodes,pt_s,subdiv_end);            %cortical principal angle threshold
disp(u); 

%% Iterate ICO-1 through ICO-End to Find Suitable Patches
tic; patchno = []; 
for subdiv = 1:min(subdiv_end,4)    
    %% Relevant Greedy Parameters
    p = eigmodes(subdiv);                                                   %read in # modes for this subdiv

    %% Form lead field matrices (highest eigenmodes of patch fields, superficial sources removed)
    [Gn,Gtheta,KK,currstr] = disjoint_leads(ico,p,pt_s,patchno,source,subdiv,subdiv_end,curr_den_fname);
    if scale_currstr == 0
        currstr = mean(currstr)*ones(size(currstr));                        %uniform scaling by mean currstr - representative strength for cortex
    end
    Gtheta = Gtheta*diag(currstr);                                          %scale cortical G fields by corresponding current strengths
            
    %% Compute Greedy Subspace Pursuit Estimates on Modes of ICO(subdiv) Patches
    disp(['Subspace pursuit in ico-' num2str(subdiv) 'p source space'])
    [~,~,trimP,Xls] = subspace_pursuit(alg,C,Gn,Gtheta,L,N,p,patchno,subdiv,stop,t1,t2,u,Y,dub,SNR,currstr);   
    X_hat{subdiv} = Xls{end};
    isequal(L, length(find(sum(X_hat{subdiv},2))))                          %check to ensure that we are finding L modes
    if isequal(subdiv,1)
        patchno = trimP;
    else
        patchno = KK(trimP);
    end

    %% Save Estimates Projected from Reduced to Hypersampled Source Space - STC and Time Series
    patch_est{subdiv} = patchno;
    tic
    plot_est{subdiv} = full_recon(alg,ico,fwd_meg,p,patchno,prefix,source,stc,subdiv,trimP,X_hat{subdiv},L,curr_den_fname,scalarX);
    toc
    for ii = 1:length(patchno)
        Xest_series = plot_est{subdiv}(ii).dip_res;                         %scalarX if dip_resovern
        Xest(subdiv).patch{ii} = Xest_series;
        Xdip(subdiv).alg = 'greedy';
    end
end
tSPIGH = toc;

%% Save and Plot Results
lk = patchno; disp(strcat('final cortical solution:  ', num2str(lk)));
savename = strcat(prefix,'-',meastype,'-',datatype,'_gr-cortest',append,'.mat');
save(savename, 'prefix','meastype','datatype','patch_est','plot_est','Xest',...
     'alg','stop','u','C','lk','subdiv','subdiv_end','-v7.3'); 
% for subcortical stage take all modes within ico(subdiv_end).patch(lk)
cat004_plot_cortical_timecourses_megeeg;
diary off;