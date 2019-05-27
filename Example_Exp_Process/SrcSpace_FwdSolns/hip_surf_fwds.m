%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%Read, Preprocess and Split Fwd Solutions for Hippocampal Surface SourceSpaces%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Read Fwd Solution
subdiv = 1;                                                                 %only ico-1 if fwd_type all, else dummy for reading src, fwd, hico

megfwd_fname = strcat(datapath,'/fwd/',prefix,'_hipsurf-meg-',fwd_type,'-fixed-fwd.fif');
eegfwd_fname = strcat(datapath,'/fwd/',prefix,'_hipsurf-eeg-',fwd_type,'-fixed-fwd.fif');
megfwd_mindist_fname = strcat(datapath,'/fwd/',prefix,'_hipsurf-meg-',fwd_type,'-mindist-fixed-fwd.fif');
eegfwd_mindist_fname = strcat(datapath,'/fwd/',prefix,'_hipsurf-eeg-',fwd_type,'-mindist-fixed-fwd.fif');
llabs = strcat(prefix,'_hipsurf-meg-',fwd_type,'-mindistout-lh-label.txt');
rlabs = strcat(prefix,'_hipsurf-meg-',fwd_type,'-mindistout-rh-label.txt');
if strcmp(meastype,'meg')
    fwd_meg = mne_read_forward_solution(megfwd_fname);                      %Hippocampal ico MEG fwd solution
    fwd_meg_mindist = mne_read_forward_solution(megfwd_mindist_fname);      %Hippocampal ico mindist fwd solution
    mindistout_labels(prefix, fwd_meg, fwd_meg_mindist,llabs, rlabs);       %Labels < mindist
    rawG = fwd_meg.sol.data;                                                %Raw G matrix to be split
    whiteG = megCw*rawG(meg_sel,:);
elseif strcmp(meastype,'eeg')
    fwd_eeg = mne_read_forward_solution(eegfwd_fname);                      %Hippocampal ico EEG fwd solution
    fwd_eeg_mindist = mne_read_forward_solution(eegfwd_mindist_fname);      %Hippocampal ico mindist fwd solution
    mindistout_labels(prefix, fwd_eeg, fwd_eeg_mindist,llabs, rlabs);       %Labels < mindist
    rawG = fwd_eeg.sol.data;                                                %Raw G matrix to be split
    whiteG = eegCw*rawG(eeg_sel-306,:);
else
    fwd_meg = mne_read_forward_solution(megfwd_fname);                      %Hippocampal ico MEG fwd solution
    fwd_eeg = mne_read_forward_solution(eegfwd_fname);                      %Hippocampal ico EEG fwd solution
    fwd_meg_mindist = mne_read_forward_solution(megfwd_mindist_fname);      %Hippocampal ico mindist fwd solution
    fwd_eeg_mindist = mne_read_forward_solution(eegfwd_mindist_fname);      %Hippocampal ico mindist fwd solution
    mindistout_labels(prefix, fwd_meg, fwd_meg_mindist, llabs, rlabs);      %Labels < mindist
    rawG = cat(1,fwd_meg.sol.data, fwd_eeg.sol.data);                       %Concatenate - Raw G matrix to be Split
    whiteG = cat(1, megCw*megP*fwd_meg.sol.data(meg_sel,:), eegCw*eegP*fwd_eeg.sol.data(eeg_sel-306,:));
    %whiteG is for calculation of whitened singular values - exclude bad channels to match dims of megCw and eegCw
end

%% Read Source Locations and Split Fwd Solution into Subdivisions 
svds_fname = strcat(datapath,'/fwd/',prefix,'_',meastype,'_hipsurf_',fwd_type,'_raw_SVDs',append,'.mat');
p_sing = 6;                                                                 %#Singular values to compute per patch 
pt_s = 10;                                                                  %Perc. superficial dipoles tolerable
nthresh = 95;                                                               %NMRA threshold for #modes per patch
if exist(svds_fname,'file') ~= 2 
  %note fwd_meg is only used for its src structures - real measurement info is in rawG
  if ~set_oct 
      cpatch_decomp(strcat(prefix,'_hipsurf'), fwd_meg, rawG, nthresh, p_sing, pt_s, subdiv, svds_fname, meastype);
  else
      cpatch_decomp_oct(strcat(prefix,'_hipsurf'), fwd_meg, rawG, nthresh, p_sing, pt_s, fwd_type, svds_fname, meastype); 
  end
end
load(svds_fname);                                                           %SVD file - used for patch decomp/structure not for whitened SVD

% Left vs. Right ICO, Source and Whitened Forward Soln
lr(1) = 'l'; lr(2) = 'r';                                                   %Left vs. right
src(1) = source(subdiv).hemisph(1);                                         %Source params for left hemisph
src(2) = source(subdiv).hemisph(2);                                         %Source params for right hemisph
puse = length(source(subdiv).hemisph(1).pinfo);                             %#patches in left hemisph
hico(1).patch = ico(subdiv).patch(1:puse);                                  %Fwd solutions for left hemisph
hico(2).patch = ico(subdiv).patch(puse+1:end);                              %Fwd solutions for right hemisph
fwd(1).rawG = rawG(:,fwd_meg.src(1).vertno);      fwd(2).rawG = rawG(:,fwd_meg.src(2).vertno);
fwd(1).whiteG = whiteG(:,fwd_meg.src(1).vertno);  fwd(2).whiteG = whiteG(:,fwd_meg.src(2).vertno);
clear fwd_meg fwd_eeg fwd_meg_mindist fwd_eeg_mindist;

%% Determine #Subdivisions and Store into Structure Format
for i = l+1 : l+2
    max_whiteS{i} = []; hipcurr_str{i} = [];  hipcurr_realstr{i} = [];      %Initialize
    % Store Relevant Information by Region
    sreg_fwd(i).reg_name = strcat(lr(i-l),'hipsurf');                       %Store the Region Name
    sreg_fwd(i).index = i;                                                  %Store the Region Index
    sreg_fwd(i).rawG = [hico(i-l).patch(:).G];                              %Full Forward Solution for this Region
    sreg_fwd(i).srclocs = [src(i-l).pinfo{:}];                              %Source Vertex Locations for this Region
    sreg_fwd(i).srcind = [hico(i-l).patch(:).lk];                           %Source Vertex Indices for this Region
    sreg_fwd(i).ndipoles = size(sreg_fwd(i).srclocs,1);                     %# Dipoles in Region
    sreg_fwd(i).volume = [];                                                %Volume of Region
    sreg_fwd(i).area = sum([hico(i-l).patch(:).A]);                         %Total Surface Area of Region 

    % Store Relevant Information by Subdivision
    num_patch = length(hico(i-l).patch);                                    %number of patches generated 
    for j = 1:length(hico(i-l).patch)
        sdiv_fwd(k).reg_name = strcat(lr(i-l),'hipsurf');                   %Store the Region Name
        all_regions{k} = sdiv_fwd(k).reg_name;                              %For Future Reference
        sdiv_fwd(k).index = j;                                              %Subdivision Index
        sdiv_fwd(k).rawG = hico(i-l).patch(j).G;                            %Subdivided Raw Fwd Soln
        sdiv_fwd(k).srclocs = round(src(i-l).rr([src(i-l).pinfo{j}]',:)*1000); %Source Vertex Locations for this Subdivision
        sdiv_fwd(k).srcind = hico(i-l).patch(j).lk;                         %Vertex #s/Source Indices (G Indices) for this Subdivision
        sdiv_fwd(k).ndipoles = size(sdiv_fwd(k).rawG,2);                    %#Dipoles in this Subdivision
        sdiv_fwd(k).volume  = [];                                           %Volume of this Subdivision
        sdiv_fwd(k).realarea = hico(i-l).patch(j).A;                        %Surface Area of this Subdivision
        sdiv_fwd(k).realcurrstr = sdiv_fwd(k).realarea*subc_cd(i);          %physical current strength
        sdiv_fwd(k).area = sreg_fwd(i).area/num_patch;                      %approx. area per patch
        if strcmp(meastype,'meg')
		whiteG = megCw*sdiv_fwd(k).rawG(meg_sel,:);                 %whiten the Split Forward Matrix
	elseif strcmp(meastype,'eeg')	
		whiteG = eegCw*sdiv_fwd(k).rawG(eeg_sel,:);	            %Whiten the Split Forward Matrix
	else
		whiteG = cat(1, megCw*sdiv_fwd(k).rawG(meg_sel,:), eegCw*sdiv_fwd(k).rawG(eeg_sel,:)); %Whiten the Split Forward Matrix
	end
        max_whiteS{i} = [max_whiteS{i} max(svd(whiteG,0))];                 %Max Singular Value = Norm
        hipcurr_str{i} = [hipcurr_str{i}, subc_cd(i)*sdiv_fwd(k).area];     %geometric current strength
        hipcurr_realstr{i} = [hipcurr_realstr{i} sdiv_fwd(k).realcurrstr];  %real geometric current strength
        k = k + 1;
    end
    
    % Compute Avg. Current Strengths Per Patch and Compare to Cortical Current Strength
    proj_currstr(i) = mean(max_whiteS{i}.*hipcurr_str{i});                  %on average current strength projected to sensor space
    norm_diff_mean(i) = proj_currstr(i)/cort_currstr_mean - 1;              %Ratio with Cortical Current Strength
    disp(strcat('subcortical current strength minus cortical strength is ...', num2str(norm_diff_mean(i)*100),'%'));
end