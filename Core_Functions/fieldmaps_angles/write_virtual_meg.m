function write_virtual_meg(evoked_fname, savename)

%% Load Files
dpath = '/autofs/eris/purdongp/users/pavitrak/matlab/mne/share/mne/mne_analyze/';
s = mne_read_bem_surfaces(strcat(dpath,'306m.fif'),true);
info = fiff_read_meas_info(evoked_fname);
evoked = fiff_read_evoked_all(evoked_fname);

%% Channel Characteristics
nChan=size(s.rr,1);                         %new # channels
info.nchan=info.nchan-2;                    %set number of channels in info structure
for iChan=1:size(info.ch_names)
    info.chs(iChan).scanno=int32(iChan);    %continuous channel naming
end
info.ch_names(305:306) = [];                %get rid of 305th and 306th MEG channel names
info.chs(305:306) = [];                     %get rid of 305th and 306th MEG channel info
info.bads = [];                             %no need of bad channels

%% Rewrite Evoked Data
evoked.evoked.nave = 1;                     %only 1 epoch so no reduction factor for field maps 
evoked.evoked.first = 0;                    %only 1 timepoint
evoked.evoked.last = 0;                     %only 1 timepoint
evoked.evoked.comment = 'Virtual MEG System';%transparent naming of file type
evoked.evoked.times = 0;                    %no times
evoked.evoked.epochs = zeros(size(evoked.evoked.epochs,1),1); %no data
evoked.evoked.epochs(305:306,:) = [];       %get rid of 305th and 306th MEG data
info.projs = struct([]);                    %ensure no projectors are loaded that will cause rotation

%% Change MEG Channel Type, Location, Orientation and Naming
for iChan=1:nChan
    info.chs(iChan).cal = double(1);        %set calibration to equal value for all channels
    
    info.chs(iChan).unit = int32(112);      %based on fiff_define_constants
    info.chs(iChan).coil_type = int32(4001);%248-Channel Whole-head System (Magnes 3600 WH)
    % http://www.lin-magdeburg.de/en/special_labs/non-invasive_brain_imaging/methods/4/index.jsp
    
    info.chs(iChan).loc(1:3,1) = s.rr(iChan,:)'; %Location
    nn = s.nn(iChan,:)';                    %Orientation
    [u, ~,~ ] = svd(eye(3,3) - nn*nn');     %from mne_add_coil_defs
    if nn'*u(:,3) < 0
        u = - u;
    end
    info.chs(iChan).loc(4:12,1) = u(:);     %Rotated Orientations
    info.chs(iChan).ch_names = ['MEG ', num2str(iChan,'%03d')]; %Sequential channel naming
end

%% Write Alias File for MNE_Rename_Channels
for i = 1:nChan
    old_chnames{i} = info.ch_names(i);
    new_chnames{i} = info.chs(i).ch_names;
    col_chnames(i,:) = strcat(old_chnames{i}, ':', new_chnames{i});
end
disp(char(col_chnames));

%% Update Evoked File with New Info Structure
data.info = info;                           %New Info File
data.evoked = evoked.evoked;                %New Evoked File
fiff_write_evoked(savename,data);           %Write New Evoked File

%% Ensure Helmet Not Shifted
figure, set(gcf,'color','white');
scatter3(s.rr(:,1),s.rr(:,2),s.rr(:,3))
hold on
for i=3:3:306
    scatter3(info.chs(i).loc(1,1),info.chs(i).loc(2,1),info.chs(i).loc(3,1),'r+');
end

end