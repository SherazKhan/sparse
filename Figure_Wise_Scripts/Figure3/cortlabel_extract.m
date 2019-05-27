%% Specify Paths, Source Space and Merged Label Files
clear; clc; close all
prefix = 'SM04'; %input('subject prefix:   ');                              %SM04, cat_004, manny
meastype = 'meg'; %input(strcat('meg_eeg','meg','eeg'));                    %meg or meg_eeg
add_allpaths;
ico_num = 2; src = mne_read_source_spaces(strcat('SM04-ico-',num2str(ico_num),'p-src.fif'));
label_fnames = {'lh.somerp1.label','rh.somerp1.label',...
    'rh.front_inf1.label','rh.front_inf2.label',...
    'lh.auditory1.label','rh.auditory1.label','lh.auditory2.label','rh.auditory2.label'};%,...
    %'lh.auditory3.label','rh.auditory3.label','lh.auditory4.label','rh.auditory4.label',...
    %'lh.frontal1.label','rh.frontal1.label',...
    %'lh.frontal2.label','rh.frontal2.label','lh.frontal3.label','rh.frontal3.label','lh.frontal4.label',...
    %'rh.frontal4.label','lh.frontal5.label','rh.frontal5.label','lh.frontal6.label','rh.frontal6.label',...
    %'lh.frontal7.label','rh.frontal7.label'};
thresh = 0.90;                                                              % Require > 90% of Patch Intersect with Chosen Regions
lh_patch_label = []; rh_patch_label = []; lh_patchno = []; rh_patchno = []; left_right = [];
circ_choice = 'som_erp';

for lln = 1:length(label_fnames)
%% Read Label
input_label = mne_read_label_file(label_fnames{lln});                       % Merged Label File
if strcmp(label_fnames{lln}(1),'l')
    lr = 1;                                                                 % Left Hemisphere
elseif strcmp(label_fnames{lln}(1),'r');                     
    lr = 0;                                                                 % Right Hemisphere
end
    
%% Initialization
if lr == 1
    lh_label_vert = input_label.vertices + 1;                               % Initialize Label Vertices
    lh_patch_ind = false(length(src(1).pinfo), 1);                          % Initialize Patch Indices
elseif lr == 0
    rh_label_vert = input_label.vertices + 1;                               % Initialize Label Vertices
    rh_patch_ind = false(length(src(1).pinfo), 1);                          % Initialize Patch Indices
end

%% Identify Patch Label Numbers and Patch Indices for Chosen Regions
for ii = 1:length(src(1).pinfo)
    if lr == 1                                                              % Left Hemisphere
        lh_int = intersect(src(1).pinfo{ii}, lh_label_vert);
        if ~isempty(lh_int) && length(lh_int)/length(src(1).pinfo{ii}) > thresh
            lh_patch_ind(ii) = true;
        end
    elseif lr == 0                                                          % Right Hemisphere
        rh_int = intersect(src(2).pinfo{ii}, rh_label_vert);
        if ~isempty(rh_int) && length(rh_int)/length(src(2).pinfo{ii}) > thresh
            rh_patch_ind(ii) = true;
        end
    end
end
if lr == 1
    lh_patch_label{lln} = src(1).vertno(lh_patch_ind) - 1;                  % Patch Label Numbers
    lh_patchno{lln} = find(lh_patch_ind);                                   % Patch Indices
elseif lr == 0
    rh_patch_label{lln} = src(2).vertno(rh_patch_ind) - 1;                  % Patch Label Numbers
    rh_patchno{lln} = find(rh_patch_ind) + length(src(1).pinfo);            % Patch Indices
end
left_right{lln} = lr;                                                       % save left right designation
    
clear input_label lh_label_vert lh_patch_ind rh_label_vert rh_patch_ind lh_int rh_int
end
save(strcat(pwd, '/', circ_choice, '_circ_results/','cregion_grps.mat'),...
'label_fnames','thresh','ico_num','left_right','lh_patch_label','lh_patchno','rh_patch_label','rh_patchno');