function read_save_sfwd(prefix, meastype, datatype, ico_num, append, fwd_type)

%% Paths and Prefixes
clc; close all; rng('default');
dpath = '/autofs/eris/purdongp/users/pavitrak/';                            %base path
mnetype = 'nightly'; disp(strcat('Using ...', mnetype, '... Version of MNE Software'));
datapath = strcat(dpath,'sourceloc_data/',datatype,'/',prefix); 
addpath(strcat(datapath, '/src/')); addpath(strcat(datapath,'/mri/'));
codepath = strcat(dpath,'sourceloc_scripts/ScSPIGH'); addpath(genpath(codepath));
add_on_name = []; %add_on_name = '_bs-full';                                %if bs full add on names to saving structure
k = 1;                                                                      %Initialize fwd solution index              

%% Regions of Interest, Current Densities and Cortical ICO to Match
sregname = {'lh','la','lp','lc','lt',...
            'rh','ra','rp','rc','rt',...
            'bsred','lic','ric',...
            'lhipsurf','rhipsurf'};                                         %'bs'
current_densities(sregname, datapath);                                      %set subcortical current density info
load(strcat(datapath,'/fwd/curr_den.mat'));                                 %save subcortical current density info
l = length(sregname)-2; 
perc_match = 0.20*(ico_num==2) + 0.05*(ico_num==3);                         %cortical ico division level and matching %age
cov_currstr_lim = 0.20;                                                     %current strength variability
if strcmp(fwd_type,'oct-3')
    set_oct = 1;                                                            %hippocampal srcspace 'oct-3'
else
    fwd_type = 'all';  set_oct = 0;                                         %hippocampal srcspace 'all' (ico-1); 
end
geom_currstr = [];                                                          %initialize vector of geometric current strengths per subdiv
diary(strcat(prefix, '_scort_ico',num2str(ico_num), '_', meastype, '_srcspace_matlog',append));

%% Cortical Areas and Current Strengths
cort_surf_fwds;                                                             %compute cortical SVDs
cpatch = ico(ico_num).patch; A_bar = mean([cpatch(:).A]);                   %mean cortical current density
for i = 1:length(cpatch)
    %csing_mean(i) = mean(diag(cpatch(i).S(1:eigmodes,1:eigmodes)));        %mean cortical singular value (after whitening)
    csing_max(i) = max(diag(cpatch(i).S(1:eigmodes,1:eigmodes)));           %max cortical singular value (after whitening)
end
cort_currstr_mean = cort_cd*A_bar*mean(csing_max);                          %mean-max cortical current strength
%cort_currstr_max = cort_cd*A_bar*max(csing_max);                           %max-max cortical current strength

%% Read, Preprocess and Split Forward Solutions for Regions Having Volume Source Spaces
for i = 1:length(sregname)-2
    if mod(i,5) == 1
        rng('default');                                                     %enforce right/left identicalness to extent possible
    end
    
    % Load Forward Solution for this Region
    reg_name = sregname{i};                                                 %Region Name
    if strcmp(meastype,'meg') || strcmp(meastype,'eeg')
        fname = strcat(datapath,'/fwd/',prefix,'_',reg_name,'_',meastype,'-mindist-fwd.fif');
        fwd = mne_read_forward_solution(fname);                             %for calculation of whitened singular values
        lfwd = size(fwd.sol.data, 2); ind1 = 1:9:lfwd;   ind2 = 2:9:lfwd;   ind3 = 3:9:lfwd; 
        ind = sort([ind1 ind2 ind3]);   clear lfwd ind1 ind2 ind3;          %Read in Only Non-Repeat Portions of Fwd Soln
        G = fwd.sol.data(:,ind);                                            %Forward Solution to Consider for Splits and Analyses
        whiteG = Cw*G(sel,:);                                               %Whitened Forward Solution to Consider for Splits and Analyses
    else
        megfname = strcat(datapath,'/fwd/',prefix,'_',reg_name,'_meg-mindist-fwd.fif');
        eegfname = strcat(datapath,'/fwd/',prefix,'_',reg_name,'_eeg-mindist-fwd.fif');
        megfwd = mne_read_forward_solution(megfname);                       %MEG Soln
        eegfwd = mne_read_forward_solution(eegfname);                       %EEG Soln
        lfwd = size(megfwd.sol.data, 2); ind1 = 1:9:lfwd;   ind2 = 2:9:lfwd;   ind3 = 3:9:lfwd; 
        ind = sort([ind1 ind2 ind3]);   clear lfwd ind1 ind2 ind3;          %Read in Only Non-Repeat Portions of Fwd Soln
        megfwd.sol.data = megfwd.sol.data(:,ind); 
        eegfwd.sol.data = eegfwd.sol.data(:,ind);
        G = cat(1,megfwd.sol.data, eegfwd.sol.data);                        %Concatenate
        whiteG = cat(1, megCw*megP*megfwd.sol.data(meg_sel,:),eegCw*eegP*eegfwd.sol.data(eeg_sel-306,:)); 
    end
    
    % Load Source Locations and Volume Information
    srclocs = importdata(strcat(datapath,'/src/',prefix,'_',reg_name,'.txt')); %Source Locations
    load(strcat(datapath,'/src/',prefix,'_',sregname{i},'_geom.mat'));      %Source volume and #dipoles                              
    
    % Check Size Consistency
    if abs(size(srclocs, 1) - size(G,2)) > 0
        disp(strcat('please run the following command in terminal for', reg_name));
        disp(['setenv reg ' reg_name]);
        disp(['mne_forward_solution --meg --accurate --src $SUBJECTS_DIR/$SUBJECT/src/$SUBJECT"_"$reg"-src.fif"',...
             ' --meas $SUBJECTS_DIR/$SUBJECT/meas/abr_plr_runall/abr_runall-',append,'-bc-ave.fif',...  %/preproc_rawdata/ABR_run1_raw.fif"
             '--mri $SUBJECTS_DIR/$SUBJECT/fwd/$SUBJECT"-ABR-bc-ave-trans.fif"',...
             ' --bem $SUBJECTS_DIR/$SUBJECT/bem/$SUBJECT"-5120-5120-5120-bem-sol.fif"' ...
             '--all --fwd $SUBJECTS_DIR/$SUBJECT/fwd/$SUBJECT"_"$reg"_meg-fwd.fif"']);
        pause;
        megfwd_nomindist = mne_read_forward_solution(strcat(datapath, '/fwd/', prefix, '_',reg_name, '_meg-fwd.fif'));
        src = mne_read_source_spaces('/autofs/eris/purdongp/users/pavitrak/sourceloc_data/',datatype,'/',prefix,'/src/',prefix,'_bs-src.fif');
        %[~, lk, ~] = intersect(src.vertno',megfwd_nomindist.src.vertno','rows');%indices of sources that are excluded in fwd_nomindist due to dist. inner skull
        [~, lk, ~] = intersect(src.vertno',megfwd.src.vertno','rows');      %G is based on fwd.source_Rr ==> select srclocs with same source locations
        %lk = megfwd_nomindist.src.vertno(logical(megfwd.src.inuse));       %ID Vertices that are in both Fwd and Mindist Fwd
        srclocs = srclocs(lk,:);                                            %doublecheck -- src2 = src.vertno(logical(fwd.src.inuse)); isequal(src2, src.vertno(lk))
    end
    clear fwd megfwd eegfwd ind 
    
    % Store Relevant Information by Region
    sreg_fwd(i).reg_name = reg_name;                                        %Store the Region Name
    sreg_fwd(i).index = i;                                                  %Store the Region Index
    sreg_fwd(i).rawG = G;                                                   %Full Forward Solution for this Region
    sreg_fwd(i).srclocs = srclocs(:,1:3);                                   %Source Locations for this Region
    sreg_fwd(i).srcind = (1:length(sreg_fwd(i).srclocs(:,1)))';             %Source Indices for this Region
    sreg_fwd(i).ndipoles = ndipoles;                                        %# Dipoles in Region
    sreg_fwd(i).volume = volume;                                            %Volume of Region
    sreg_fwd(i).area = [];                                                  %Surface Area of Region
    
    % Number of Subdivisions per Regions Having Volume Source Spaces
    [Gw, Uw, Sw, Vw, pw, nmw, ~] = nmra(whiteG, 0, 1, 6, 95);               %compute singular values of whitened full region fwd
    ssing(i) = max(diag(Sw(1:pw,1:pw)));                                    %max singular value
    sreg_fwd(i).whiteS = Sw; sreg_fwd(i).numwhite_modes = pw;               %relevant SVD characteristics
    V(i) = sreg_fwd(i).volume;                                              %volume of this region
    scort_currstr(i) = subc_cd(i)*V(i)*ssing(i);                            %current strength of this region
    if ~strcmp(reg_name,'lic') && ~strcmp(reg_name,'ric')
    num_divs(i) = (scort_currstr(i)/cort_currstr_mean)^(2/3);               %number of subdivisions of this region
    disp(scort_currstr(i)); disp(num_divs(i));
    if isempty(geom_currstr)
        geom_currstr = subc_cd(i)*V(i)/num_divs(i);                         %benchmark geometric current strength
    end
    [ssing(i),  scort_currstr(i), num_divs(i), split_ind, split_srclocs, geom_currstr] = ...
    check_numdivs(srclocs, num_divs(i), strcat(datapath, '/src/',reg_name), 0, G, ...
    cort_currstr_mean, subc_cd(i), V(i), ndipoles, geom_currstr, perc_match, cov_currstr_lim, megCw*megP, meg_sel, eegCw*eegP, eeg_sel);
    else
    num_divs(i) = 1; ssing(i) = 0; scort_currstr(i) = 0; 
    [split_ind, split_srclocs] = split_subdiv(srclocs, num_divs(i), strcat(datapath, '/src/',reg_name), 0);
    end
    
    % Split Fwd Solution into Subdivisions and Store Relevant Information by Region  
    for j = 1:length(split_ind)
        sdiv_fwd(k).reg_name = reg_name;                                    %Store the Region Name
        all_regions{k} = sdiv_fwd(k).reg_name;                              %For Future Reference
        sdiv_fwd(k).index = j;                                              %Subdivision Index
        sdiv_fwd(k).rawG = G(:,split_ind{j});                               %Subdivided Raw Fwd Soln
        sdiv_fwd(k).srclocs = split_srclocs{j}(:,1:3);                      %Source Locations for this Subdivision
        sdiv_fwd(k).srcind = find(split_ind{j});                            %Source Indices (G Indices) for this Subdivision
        sdiv_fwd(k).ndipoles = floor(size(sdiv_fwd(k).rawG,2)/3);           %#Dipoles in this Subdivision
        sdiv_fwd(k).realvol = volume/ndipoles*floor(size(split_ind{j},2)/3);%Physical volume
        sdiv_fwd(k).realcurrstr = sdiv_fwd(k).realvol*subc_cd(i);           %Physical current strength
        sdiv_fwd(k).volume  = V(i)/num_divs(i);                             %Average Volume of this Subdivision
        sdiv_fwd(k).area = [];                                              %Surface Area of this Subdivision
        k = k + 1;
    end
end

%% Read, Preprocess and Split Fwd Solutions for Hippocampal Surface SourceSpaces
clear reg_name split_ind split_srclocs G ind fwd
hip_surf_fwds;
scort_currstr(i-1) = unique(hipcurr_str{i-1});	scort_currstr(i) = unique(hipcurr_str{i});

%% Save Fwds for Future Access
% sdiv_fwd and sreg_fwd have reg_name/index, rawG, src locations, volume, #dipole fields
save(strcat(datapath,'/fwd/',prefix,'_subc_',meastype,'_plr-prewhite_volumes',add_on_name,append,'.mat'),...
     'sregname', 'sdiv_fwd','sreg_fwd','ssing','csing_max','num_divs',...
     'scort_currstr','cort_currstr_mean','ico_num','A_bar');
disp('completed saving subdivided source space');
disp(num_divs);
nn_diff = num2str(abs(scort_currstr - cort_currstr_mean)/cort_currstr_mean*100); disp(nn_diff);
diary off;

end