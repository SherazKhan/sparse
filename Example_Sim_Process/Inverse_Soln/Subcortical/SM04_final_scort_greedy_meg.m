%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SPIGH SOLUTION TO THE NEUROMAGNETIC INVERSE PROBLEM
%Localize Subcortical Sources with Subspace Pursuit (Reduced State-Space Model)
%Written by Pavitra Krishnaswamy (pavitrak@nmr.mgh.harvard.edu)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Set Paths, Load Cortical Forward Solutions, Meas Files
clear; clc; close all; 
prefix = 'SM04'; meastype = 'meg'; datatype = 'SEP'; append = '_overlap';
add_allpaths;                                                               %add all paths
diary(strcat(prefix, '_', meastype, '_scort_matlog',append));
s = strfind(pwd,'/'); cortpath = pwd; cortpath = strcat(cortpath(1:s(end)), '160720_final_SM04_cortsoln_ico3/');
cortest_fname = strcat(cortpath, prefix,'-',meastype,'-',datatype,'_gr-cortest',append,'.mat');%_L8, _SNR10
csvds_fname = strcat(datapath, 'fwd/', prefix,'_cort_',meastype,'_prewhite','_SVDs',append,'.mat'); 
curr_den_fname = strcat(datapath, 'fwd/curr_den.mat');                      %for cortical current density
meas_fname = strcat(datapath, 'meas/', prefix,'_vpl_somato_meg_sim',append,'.mat'); %measurement file with all info
subc_svds_fname = strcat(datapath, 'fwd/', prefix,'_subc_',meastype,'_prewhite','_SVDs',append,'.mat'); %subcortical SVD file
saveon = 1; 

%% Cortical Soln: Load Files and Parameter Settings
load(cortest_fname,'alg','patch_est','stop','subdiv_end','u');              %Load cortical estimates
cdiv = subdiv_end;                                                          %What cortical subdiv to consider for hybrid source space
load(csvds_fname,'ico','p'); eigmodes = p;  clear p;                        %forward solutions: patch decompositions
lk = patch_est{cdiv};                                                       %chosen cortical patches for this cdiv
ico = ico(cdiv);                                                            %the ico we need to consider
cp = eigmodes(cdiv);                                                        %# eigenmodes for cortical stage
load(curr_den_fname,'cort_cd');                                             %[nAm/mm^2] - cortical current density
clkmax = length(ico.patch); lkmax = clkmax;                                 %this is for the cortical subdiv we consider here

%% Subcortical Soln: Load Whitened Measurement, Covariance and Compute Principal Angle Threshold
load(meas_fname,'Cw','Yw','C','N','t1','t2','stc','cij','scalarX','time'); %'Fs' %meas, covariance, whitener for subcortical step
Y = Yw;                                                                     %use whitened data
incl_sregions = {'la','lp','lc','lt','ra','rp','rc','rt','lhipsurf','rhipsurf','bsred'}; %all regions to include in the inverse solution
pa = scort_mutcoh_neighbors(subc_svds_fname); %0.95; %0.88;                 %average worst case mutcoh b/w subdivisions in same anatomic region - deep
disp(['subcortical neighborhood mutual coherence...', num2str(pa)]);        %display principal angle threshold 

%% Hybrid Solution: Set Parameters, Time Points of Interest and Offset
fullcort = 0;                                                               %0 hierarchial greedy, 1 for no hierarchy just sparsity
snummodes = 0;                                                              %0 for all modes, 1 for top mode
select_time = 0;                                                            %select full time range for evoked response
if select_time == 1                                                         %select known times for evoked response
    t1 = 1; t2 = 45;                                                        %SM04 - 1st 15 msec of ERP
end
offset = 0;                                                                 %inconsequential param for time vector
dub = 0;                                                                    %inconsequential param for SP run: doubling heuristic
if strcmp(append, '_overlap')
    L = 12;                                                                 %target sparsity level
else    
    L = 13;%14;                                                             %target sparsity level
end
sSNR = 16; %10;                                                             %typical ERP case power SNR 9, can also choose 16

%% Cortical Fwd Solution: Lead Fields for Best Cortical Patches
disp('Obtaining Forward Solutions for Best Cortical Regions');
Gcn = []; Gctheta = []; Sc = []; A = []; clk = [];
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
lk = clk;                                                                   %cortical patches to go into subcortical soln

%% Subcortical Fwd Solution: Load and Concatenate with Cortical Fwds
disp('Adding in Forward Solutions for Subcortical Regions');
%Gn = []; Gtheta = []; lk = [];                                             %to run purely subcortical soln w/o any cortical components
load(subc_svds_fname,'sdiv_fwd','sregname');                                %subcortical fwds by subdivisions with regnames
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
    
    %Load in Fwd Solutions of ith Region
    Gstheta = [sfwd(ind_sfwd).Gstheta_p];                                   %all the Gstheta for ith region
    Gsn = [sfwd(ind_sfwd).Us_p];                                            %all the Us for ith region
    Gtheta = cat(2,Gtheta, Gstheta);                                        %put Gstheta into overall Gtheta
    Gn = cat(2, Gn, Gsn);                                                   %put Gn into overall Gn
    currstr = [currstr,sfwd(ind_sfwd).currstr];                             %put sfwd.currstr into overall currstr
    eff_currstr = [eff_currstr, sfwd(ind_sfwd).eff_currstr];                %put sfwd.efff_currstr into overall eff_currstr
    
    %Define Patch Indices for ith Region
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
Y = Y(ssel,t1:t2);                                                          %use selected part of measurement
N = size(Y,1);                                                              %#channels being used
t = size(Y,2);                                                              %should be same as t2-t1+1
time_est = offset+linspace(t1-1,t2,t)*stc(1).tstep;                         %time for plotting estimates
Gn = Gn(ssel,:);                                                            %selected rows of forward modes
Gtheta = Gtheta(ssel,:);                                                    %selected rows of forward solution
C = C(ssel,ssel);                                                           %selected rows of covariance matrix

%% Find Subcortical Estimates - Greedy and MNE
disp('Estimation of Source Currents in Joint Cortical-Subcortical SourceSpace');
% Principal Angles
num_cort = length(unique(clk));     num_scort = length(sdiv_fwd);           %no. of patches and subdivisions in hybrid source space
pa = round(sum([pa*num_scort, u(cdiv)*num_cort])/(num_scort+num_cort),2);   %for hybrid stage include cortical stages in the neighborhood
subdiv = 5; u_scort(subdiv) = pa;                                           %inconsequential param for SP run: each mode is independent for subspace pursuit

% Subspace Pursuit
%currstr = eff_currstr; Gtheta = Gn;                                        %if only want to consider patterns
[~,~,trimP,Xgreedy] = subspace_pursuit_scort(alg,C,Gn,Gtheta,L,N,lk,subdiv,stop,t1,t2,u_scort,Y,dub,sSNR, currstr);
disp(strcat('final hybrid solution:  ' , num2str(lk(trimP))));              %lk indices found by greedy

% MNE
[Xmne,~] = source_estimates(alg,C,Gtheta,N,Y,t1,t2,sSNR,stop, currstr);     %benchmark against straight MNE
proxy = Xmne; n = sqrt(dot(proxy,proxy,2)); figure, plot(n);                %MNE proxy shape for evaluation
%figure, plot(currstr); figure, plot(eff_currstr);  

%% Save Results
if fullcort == 0
    savename = strcat(prefix,'-',meastype,'-',datatype,'_scortest',append,'.mat');
else
    savename = strcat(prefix,'-',meastype,'-',datatype,'_1stage_scortest',append,'.mat');
end
if saveon
    save(savename,'cdiv','all_regions','sregnums','regname','svolume',...
    'select_time','fullcort','snummodes','pa','u_scort','dub','L','sSNR',...
    'alg', 'N', 't1', 't2', 'time_est',...
    'cp','c_curr_str','sp','currstr','ssubdiv','Gn', 'Gtheta',...
    'lk', 'trimP','Xgreedy','Xmne','clkmax','lkmax');
end
diary off;