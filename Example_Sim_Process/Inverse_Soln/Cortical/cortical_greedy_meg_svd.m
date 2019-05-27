%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SPIGH SOLUTION TO THE NEUROMAGNETIC INVERSE PROBLEM
%Localize Cortical Sources with Subspace Pursuit (Reduced State-Space Model)
%Written by Gabriel Obregon (obregon@nmr.mgh.harvard.edu)
%Modified and Adapted by Pavitra Krishnaswamy (pavitrak@nmr.mgh.harvard.edu)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Set Paths, Load Cortical Forward Solutions, Meas Files
clear; clc; close all; 
prefix = 'SM04'; meastype = 'meg'; datatype = 'SEP'; append = '_overlap'; saveon = 0;
if saveon
    diary(strcat(prefix, '_', meastype, '_cort_matlog',append));
end
datapath = strcat('/autofs/eris/purdongp/users/pavitrak/sourceloc_data/',datatype,'/', prefix,'/');
addpath(genpath('/autofs/eris/purdongp/users/pavitrak/sourceloc_scripts/ScSPIGH/')); 
cfwd_fname = strcat(datapath, '/fwd/', prefix,'-', meastype,'-all-fixed-fwd.fif');%cortical forward solution
fwd = mne_read_forward_solution(cfwd_fname);                                %read fwd soln
csvds_fname = strcat(datapath, 'fwd/', prefix,'_cort_',meastype,'_prewhite','_SVDs',append,'.mat'); 
load(csvds_fname); eigmodes = p;  clear p;                                  %forward solutions: patch decompositions
curr_den_fname = strcat(datapath, 'fwd/curr_den.mat');                      %for cortical current density
load(strcat(datapath, 'meas/', prefix,'_vpl_somato_meg_sim',append,'.mat'),...%meas, covariance, whitener
           'Cw','Yw','C','N','t1','t2','stc', 'scalarX','time','sim_on');
Y = Yw;                                                                     %use whitened data

%% Set Greedy Parameters
pt_s = 25;                                                                  %perc. superficial dipoles tolerable
alg = 'MNE'; %'dSPM'; %'sMAP-EM'                                            % if sMAP-EM do gradiometers only
if strcmp(alg,'sMAP-EM');
    stop = 7;
else
    stop = [];
end
subdiv_end = 3;                                                             % SPIGH to ICO-3 (could also run MNE to ICO-4 after)
if strcmp(append,'_overlap')
    SNR = 16; %10;                                                           % typical ERP case 3
    L = 8; dub = 0;                                                         % target sparsity level and doubling heuristic
else
    SNR = 16;                                                               % typical ERP case 3  
    L = 7; dub = 0; %L = 8;                                                 % target sparsity level and doubling heuristic
end
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
    plot_est{subdiv} = full_recon(alg,ico,fwd,p,patchno,prefix,source,stc,subdiv,trimP,X_hat{subdiv},L,curr_den_fname,scalarX);
    toc
    for ii = 1:length(patchno)
        Xest_series = plot_est{subdiv}(ii).dip_res;                         %scalarX if dip_resovern
        Xest(subdiv).patch{ii} = Xest_series;
        Xdip(subdiv).alg = 'greedy';
    end
end
tSPIGH = toc;

%% Save Results
if saveon
lk = patchno; disp(strcat('final cortical solution:  ', num2str(lk)));
savename = strcat(prefix,'-',meastype,'-',datatype,'_gr-cortest',append,'.mat');%_SNR10
save(savename, 'prefix','meastype','datatype','patch_est','plot_est','Xest',...
     'alg','stop','u','C','lk','subdiv','subdiv_end','-v7.3'); 
% for subcortical stage take all modes within ico(subdiv_end).patch(lk)
SM04_plot_cortical_timecourses_meg;
diary off;
end