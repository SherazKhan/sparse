function patch_areas(prefix, svds_fname, fwd, subdiv_stop)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Find Areas of Cortical Patches
% Written by Pavitra Krishnaswamy (pavitrak@nmr.mgh.harvard.edu)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Load SVDs File to Append Patch Areas to
load(svds_fname);
disp('patch areas computing')

for subdiv = 1:min(subdiv_stop,4)
    %% Read Forward and Rad Files
    lrad = mne_read_w_file(strcat(prefix,'-ico-',num2str(subdiv),'p-rad-lh.w'));
    lrad.data = lrad.data(logical(fwd.src(1).inuse));
    rrad = mne_read_w_file(strcat(prefix,'-ico-',num2str(subdiv), 'p-rad-rh.w'));
    rrad.data = rrad.data(logical(fwd.src(2).inuse));
    nuse = source(subdiv).hemisph(1).nuse;
    
    %% Radii for our Patches
    for ii = 1:nuse
        prad(ii) = setdiff(unique(lrad.data(ico(subdiv).patch(ii).lk),'stable'),15,'stable');
        ico(subdiv).patch(ii).A = pi*prad(ii).^2;
        prad(ii+nuse) = setdiff(unique(rrad.data(ico(subdiv).patch(ii+nuse).lk-double(fwd.src(1).nuse)),'stable'),15,'stable');
        ico(subdiv).patch(ii+nuse).A = pi*prad(ii+nuse).^2;
    end    
    
    %% Plot Areas of Patches
    figure(30), set(gcf,'color','white');
    if subdiv_stop >= 4
        subplot(2,2,subdiv), mareas = pi*prad.^2; plot(mareas,'k');
    else
        mareas = pi*prad.^2; plot(mareas, 'k');
    end
    title(strcat('Areas for Ico-', num2str(subdiv),'Patches'));
    xlabel('Patch Numbers'); ylabel('Area (sq.mm)');
end

%% Save
save(svds_fname,'ico','source','-append','-v7.3');
disp('patch areas done');
end