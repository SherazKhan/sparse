function excl_medial_labels(subjdir, prefix, ico_num, csvd_fname, thresh)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Function to exclude medial cortical labels
%Written by: Pavitra Krishnaswamy
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Set up Paths
src_path = strcat(subjdir,'/src/');                                         % for source space files
label_path = strcat(subjdir,'/label/');                                     % for labels
fwd_path = strcat(subjdir, '/fwd/');                                        % for forward solutions
src = mne_read_source_spaces(strcat(src_path, prefix, '-ico-',num2str(ico_num),'p-src.fif'));

%% Labels to Exclude
atlas_label_path = strcat(label_path, prefix, '-Destrieux/');               % If Destrieux 2005 Atlas: 2005_labels Instead
label_fnames{1} = strcat(atlas_label_path, 'lh.Unknown.label');             % If Destrieux 2005 Atlas: : lh.Medial_wall.label
label_fnames{2} = strcat(atlas_label_path, 'rh.Unknown.label');             % If Destrieux 2005 Atlas: : rh.Medial_wall.label  
if exist(label_fnames{1},'file') ~= 2 
disp('please run below commands in terminal');
disp('mri_annotation2label --subject $SUBJECT --hemi rh --annotation $SUBJECTS_DIR/$SUBJECT/label/rh.aparc.a2009s.annot --outdir $SUBJECTS_DIR/$SUBJECT/label/$SUBJECT"-Destrieux"');
disp('mri_annotation2label --subject $SUBJECT --hemi lh --annotation $SUBJECTS_DIR/$SUBJECT/label/lh.aparc.a2009s.annot --outdir $SUBJECTS_DIR/$SUBJECT/label/$SUBJECT"-Destrieux"');
pause; 
end

for lln = 1:length(label_fnames)
%% Read in Label
input_label = mne_read_label_file(label_fnames{lln});                       % Read Label File
if strcmp(label_fnames{lln}(length(atlas_label_path)+1),'l')
    lr = 1;                                                                 % Left Hemisphere
elseif strcmp(label_fnames{lln}(length(atlas_label_path)+1),'r');                     
    lr = 0;                                                                 % Right Hemisphere
end

%% Initialization
if lr == 1
    lh_label_vert = input_label.vertices + 1;                               % Initialize Label Vertices
    lh_patch_ind = false(length(src(1).pinfo), 1);                          % Initialize Patch Indices
elseif lr == 0
    rh_label_vert = input_label.vertices + 1;                               % Initialize Label Vertices
    rh_patch_ind = false(length(src(2).pinfo), 1);                          % Initialize Patch Indices
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

%% Save Results
excl_patch = double([lh_patchno{1}' rh_patchno{2}']);
save(strcat(fwd_path, 'exclude_patchlbl_ico',num2str(ico_num),'.mat'),...
    'label_fnames','thresh','ico_num','left_right',...
    'excl_patch','lh_patch_label','lh_patchno','rh_patch_label','rh_patchno');

%% Flag Patches for Exclusion in Patch Decomp Files
% Initialize
load(csvd_fname, 'ico'); 
for ii = 1:length(ico(ico_num).patch)
	ico(ico_num).patch(ii).exclude = 0;
end
% Left Hemisphere
for i = 1:length(lh_patchno{1})
	ico(ico_num).patch(lh_patchno{1}(i)).exclude = 1;
end
% Right Hemisphere
for i = 1:length(rh_patchno{2})
	ico(ico_num).patch(rh_patchno{2}(i)).exclude = 1;
end
% Resave by Appending and Clear
save(csvd_fname, 'ico', '-append');
clear ico; 

end