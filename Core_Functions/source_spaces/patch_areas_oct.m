function patch_areas_oct(prefix, svds_fname, fwd, fwd_type)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Find Areas of Cortical Patches
% Written by Pavitra Krishnaswamy (pavitrak@nmr.mgh.harvard.edu)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Load SVDs File to Append Patch Areas to
load(svds_fname);
disp('patch areas computing')

%% Read Forward and Rad Files
lrad = mne_read_w_file(strcat(prefix,'-',fwd_type,'p-rad-lh.w'));
lrad.data = lrad.data(logical(fwd.src(1).inuse));
rrad = mne_read_w_file(strcat(prefix,'-',fwd_type,'p-rad-rh.w'));
rrad.data = rrad.data(logical(fwd.src(2).inuse));
nuse = source.hemisph(1).nuse;
    
%% Radii for our Patches
for ii = 1:nuse
	prad(ii) = setdiff(unique(lrad.data(ico.patch(ii).lk),'stable'),15,'stable');
    ico.patch(ii).A = pi*prad(ii).^2;
    prad(ii+nuse) = setdiff(unique(rrad.data(ico.patch(ii+nuse).lk-double(fwd.src(1).nuse)),'stable'),15,'stable');
    ico.patch(ii+nuse).A = pi*prad(ii+nuse).^2;
end    
    
%% Plot Areas of Patches
figure(30), set(gcf,'color','white');
mareas = pi*prad.^2; plot(mareas,'k');
title(['Areas for ', fwd_type, ' Patches']);
xlabel('Patch Numbers'); ylabel('Area (sq.mm)');

%% Save
save(svds_fname,'ico','source','-append','-v7.3');
disp('patch areas done');
end