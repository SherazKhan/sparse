%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SPIGH SOLUTION TO THE NEUROMAGNETIC INVERSE PROBLEM
%Localize Subcortical Sources with Subspace Pursuit (Reduced State-Space Model)
%Written by Pavitra Krishnaswamy (pavitrak@nmr.mgh.harvard.edu)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Set Paths, Load Cortical Forward Solutions, Meas Files
clear; clc; close all; 
prefix = 'SM04'; meastype = 'meg'; datatype = 'Illustrations';
dpath = '/autofs/eris/purdongp/users/pavitrak/';  add_allpaths;             %add all paths
diary(strcat(prefix, '_', meastype, '_scort_matlog'));
savepath = strcat(pwd,'/2016_Results_MatsFigs/'); 
sinv_fname = strcat(savepath, prefix,'-',meastype,'-',datatype,'_sinv.mat');
load(sinv_fname,'Cw','alg','stop','cdiv','patch_est','cij','SNR','csvds_fname','ssvds_fname'); %ico, #divs, selected cortical patches,sdivs, greedy settings
curr_den_fname = strcat(datapath, '/fwd/curr_den.mat');                     %for cortical current density
sim_fname = strcat(savepath, prefix,'-', meastype, '-',datatype,'_noiseless_sim.mat'); %measurement file with all info
saveon = 1; ploton = 1; greedy = 1;                                         %0/1 for MNE vs. greedy for emp_resmat

%% Cortical Soln: Load Files and Parameter Settings
load(csvds_fname,'ico','p'); ico = ico(cdiv); eigmodes = p(cdiv);  clear p; %forward solutions: patch decompositions
lk = unique(patch_est);                                                     %chosen cortical patches for this cdiv
load(curr_den_fname,'cort_cd');                                             %[nAm/mm^2] - cortical current density
clkmax = length(ico.patch); lkmax = clkmax;                                 %this is for the cortical subdiv we consider here

%% Subcortical Soln: Load Whitened Measurement, Covariance and Compute Principal Angle Threshold
load(ssvds_fname,'sdiv_fwd','sregname');                                    %subcortical SVDs
load(sim_fname,'xc','xs','Y','N','t1','t2','t','C');                        %whitened meas, whitened covariance, whitener for subcortical step
incl_sregions = {'lhipsurf','la','lp','lc','lt','bsred','rt','rc','rp','ra','rhipsurf'}; %all regions to include in the inverse solution
pa = scort_mutcoh_neighbors(ssvds_fname);                                   %average worst case mutcoh b/w subdivisions in same anatomic region - deep
disp(['subcortical neighborhood mutual coherence...', num2str(pa)]);        %display principal angle threshold 

%% Hybrid Solution: Set Parameters, Time Points of Interest and Offset
fullcort = 0;                                                               %0 hierarchial greedy, 1 for no hierarchy just sparsity
snummodes = 1; cnummodes = snummodes;                                       %0 for all modes, 1 for top mode
offset = 0;                                                                 %inconsequential param for time vector
L = 1; dub = 0;                                                             %overall target spar
sSNR = SNR;                                                                 %typical ERP case power SNR 9, can also choose 16

%% Cortical Fwd Solution: Lead Fields for Best Cortical Patches
disp('Obtaining Forward Solutions for Best Cortical Regions');
Gcn = []; Gctheta = []; Sc = []; A = []; clk = [];
if cnummodes
    cp = 1;                                                                 %represent each division by top white mode 
else
    cp = eigmodes; 
end
if fullcort == 1                                                            %if inputting full cortical solution in 1 stage
    lk = 1:1:length(ico.patch);                                             %full iconum source space
end
for ii = 1:length(lk)                                                       %for each cortical patch
    Gcn = cat(2,Gcn,ico.patch(lk(ii)).U(:,1:cp));                           %modes
    Gctheta = cat(2,Gctheta,ico.patch(lk(ii)).U(:,1:cp)*ico.patch(lk(ii)).S(1:cp,1:cp)); %forward soln
    Sc = cat(2,Sc,diag(ico.patch(lk(ii)).S(1:cp,1:cp))');                   %singular values
    %A = [A repmat(ico.patch(lk(ii)).A,1,cp)];                              %individual patch areas
    clk = [clk repmat(lk(ii),1,cp)];                                        %as each mode is independent, repeat patchno for each mode
end
A = repmat(mean([ico.patch(:).A]),1,length(clk));                           %average area of any cortical patch
c_curr_str = A*cort_cd;                                                     %average current strength of a cortical patch
c_eff_currstr = c_curr_str.*Sc;
Gn = Gcn;                                                                   %U's dont scale by current strength
Gtheta = Gctheta*diag(c_curr_str);                                          %scale Gctheta by cortical current strength
currstr = c_curr_str;                                                       %current strength
eff_currstr = c_eff_currstr;                                                %effective current strength
left_right = lk < length(ico.patch)/2;                                      %if left this is 1
reg_names = cell(1, length(lk)); 
left = repmat({'lcortical'},1,length(lk));    right = repmat({'rcortical'},1,length(lk));
reg_names(left_right) = left(left_right); reg_names(~left_right) = right(~left_right); %initialize region names
lk = clk;                                                                   %cortical patches to go into subcortical soln
clkmax = length(ico.patch); lkmax = clkmax;                                 %all cortical patches have same #modes

%% Subcortical Fwd Solution: Load and Concatenate with Cortical Fwds
disp('Adding in Forward Solutions for Subcortical Regions');
%Gn = []; Gtheta = []; lk = [];                                             %to run purely subcortical soln w/o any cortical components
load(ssvds_fname,'sdiv_fwd','sregname');                                    %subcortical fwds by subdivisions with regnames
all_sregnums = 1:length(sregname);                                          %all region indices
sregnums = all_sregnums(ismember(sregname, incl_sregions));                 %region indices corresponding to included scort regions

% Read in Subcortical Fwd Structure For Easy Indexing
for k = 1:length(sdiv_fwd)
    pw(k) = sdiv_fwd(k).numwhite_modes;             
    if snummodes
    sp(k) = 1;                                                              %#white modes - max([sdiv_fwd(:).numwhite_modes])
    else
    sp(k) = pw(k); %best_white_sdiv_modes;                                  %#white modes - max([sdiv_fwd(:).numwhite_modes])
    end
   	all_regions{k} = sdiv_fwd(k).reg_name;                                  %names of regions
    sfwd(k).currstr = repmat(sdiv_fwd(k).currstr,1,sp(k));                  %current strength
    sfwd(k).Gstheta_p = sdiv_fwd(k).whiteGstheta(:,1:sp(k));                %whitened G cols (scaled by currstr)
    sfwd(k).Us_p = sdiv_fwd(k).whiteUs(:,1:sp(k));                          %whitened modes (not scaled by currstr)
    sfwd(k).whiteS = repmat(sdiv_fwd(k).whiteS(1:sp(k),1:sp(k)), 1, sp(k)); %whitened singular values - norms of whitened G cols when sp is 1
    sfwd(k).eff_currstr = sfwd(k).currstr*sfwd(k).whiteS;                   %effective current strength
end

for i = 1:length(sregnums)
    %Identify ith Region
    regname{i} = sregname{sregnums(i)};                                     %name of ith region selected
    ind_sfwd = strncmp(regname{i}, all_regions, 4);                         %find indices of sdiv_fwd that match regname{i}
    ind = find(ind_sfwd); 
    
    %Load in Fwd Solutions of ith Region
    Gstheta = [sfwd(ind_sfwd).Gstheta_p];                                   %all the Gstheta for ith region
    Gsn = [sfwd(ind_sfwd).Us_p];                                            %all the Us for ith region
    Gtheta = cat(2,Gtheta, Gstheta);                                        %put Gstheta into overall Gtheta
    Gn = cat(2, Gn, Gsn);                                                   %put Gn into overall Gn
    currstr = [currstr,sfwd(ind_sfwd).currstr];                             %put sfwd.currstr into overall currstr
    eff_currstr = [eff_currstr, sfwd(ind_sfwd).eff_currstr];                %put sfwd.efff_currstr into overall eff_currstr
    reg_names = [reg_names repmat(all_regions(ind(1)), 1, length(ind))];

    %Define Division Indices for ith Region
    numeigen(i) = size(Gstheta,2);                                          %number of eigenmodes for ith region
    slk = lkmax+1:1:lkmax+numeigen(i);                                      %lk indices for ith region
    ssubdiv{i} = find(ind_sfwd);                                            %subdivisions of region i
    nummodes = sp(ind_sfwd);                                                %#modes to use in region i's subdivisions
    c = 1;                                                                  %counter
    for j = 1:length(ssubdiv{i})                                            %length(ssubdiv{i}) is #subdivisions
        spj = nummodes(j);                                                  %#modes for ith region jth subdivision
        svolume{i,j} = slk(c:c+spj-1);                                      %indices for ith region jth subdivision
        c = c + spj;
    end
    
    %Add to Global Pool, Index Accordingly
    lk = [lk slk]; lkmax = max(lk);                                         %grow lk from which to choose 
end
disp(lk);                                                                   %preview lk of choice

%% Measurement, Forward Solution and Covariance
disp('Formatting Measurement to Study');
ssel = 1:1:N;                                                               %use both gradiometers and magnetometers
Y = Y(ssel,:);                                                              %use selected part of measurement
N = size(Y,1);                                                              %#channels being used
t = size(Y,2);                                                              %should be same as t2-t1+1
time_est = 1;                                                               %time for plotting estimates
Gn = Gn(ssel,:);                                                            %selected rows of forward modes
Gtheta = Gtheta(ssel,:);                                                    %selected rows of forward solution
C = C(ssel,ssel);                                                           %selected rows of covariance matrix

%% Find Subcortical Estimates - Greedy and MNE
disp('Estimation of Source Currents in Joint Cortical-Subcortical SourceSpace');
% Principal Angles
num_cort = length(unique(clk));     num_scort = length(sdiv_fwd);           %no. of patches and subdivisions in hybrid source space
%pa = round(sum([pa*num_scort, u(cdiv)*num_cort])/(num_scort+num_cort),2);   %for hybrid stage include cortical stages in the neighborhood
subdiv = 5; u_scort(subdiv) = pa;                                           %inconsequential param for SP run: each mode is independent for subspace pursuit

for src_ind = 1:size(Y,2)
    %currstr = eff_currstr; Gtheta = Gn;                                    %if only want to consider patterns
    
    % MNE    
    [Xmne{src_ind},~] = source_estimates(alg,C,Gtheta,N,Y(:,src_ind),t1,t2,SNR,stop,currstr);  %benchmark against straight MNE
    proxy{src_ind} = Xmne{src_ind};                                         %proxy for greedy
    norm_proxy{src_ind} = sqrt(dot(proxy{src_ind},proxy{src_ind},2)); %figure, plot(n); %MNE proxy shape for evaluation
    %figure, plot(currstr); figure, plot(eff_currstr);  
    [~, trimP_mne{src_ind}] = max(norm_proxy{src_ind});
    
    % Subspace Pursuit
    [~,~,trimP{src_ind},Xgreedy{src_ind}] = ...
    subspace_pursuit_scort(alg,C,Gn,Gtheta,L,N,lk,subdiv,stop,t1,t2,u_scort,Y(:,src_ind),dub,SNR,currstr);
    %disp(strcat('final hybrid solution:  ' , num2str(lk(trimP))));              %lk indices found by greedy
    subcort_ind{src_ind} = lk(trimP{src_ind}) > clkmax;                     %subcortical indices
    close all;
end

%% Compute Empirical Resolution Matrix
if greedy == 0
    trimP_best = trimP_mne; X_best = Xmne;                                  %sparse_cort_all_cort
else
    trimP_best = trimP;     X_best = Xgreedy;                               %greedy
end
K = emp_resmat(xc, xs, lk, clk, cp, regname, ssubdiv, svolume, Gtheta, trimP_best, X_best, ploton); 

%% Save Results
if fullcort == 0
    choose_cort = 'sparse';    
    if greedy == 0
        choose_scort = 'all_emp';                                           %mne results, not analytically computed
    else
        choose_scort = 'sparse';                                            %greedy results
    end
    savename = strcat(savepath, prefix, '_', choose_cort, '_cort_', choose_scort, '_scort_resmat.mat');
else
    choose_cort = 'all';       choose_scort = 'all';                        % sparse run on all cortical and all subcortical
    savename = strcat(savepath, prefix, '_', choose_cort, '_cort_', choose_scort, '_scort_greedy.mat');
end
if saveon
    save(savename,...
    'clkmax','lk','clk','lkmax','cdiv','sregnums','regname','all_regions',...
    'nummodes','L','dub','sSNR','pa','u_scort','alg','fullcort',...
    'reg_names','ssubdiv','slk','currstr','Gn', 'Gtheta','svolume',...
    'N', 't1', 't2', 'Y', 'time_est',...
    'trimP','Xgreedy','trimP_mne','Xmne','norm_proxy','proxy','subcort_ind','K');
end
diary off;