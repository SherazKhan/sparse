% Mesh wtih Patch Locatons
% Written by Gabriel Obregon-Henao (obregon@nmr.mgh.harvard.edu)
% Edited by Pavitra Krishnaswamy (pavitrak@nmr.mgh.harvard.edu)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Paths and Parameters 
close all; clear all; clc;
prefix = 'SM04'; meastype = 'meg'; 
addpath(genpath(strcat('/autofs/cluster/purdonlab/users/pavitrak/sourceloc_data/Illustrations/',prefix))); 
fwd = mne_read_forward_solution(strcat(prefix,'-',meastype,'-all-fixed-fwd.fif'));
src = mne_read_source_spaces(strcat(prefix,'-ico-2p-src.fif'));
ico_2_labelno = [58397, 80245, 28238, 84729]; %cs1, cs2, ppc, iS2
nuse = 162;
n2 = [];

%% ICO-2 Patches
for i = 1:length(ico_2_labelno)
    % Read in Selected Patch
    [~,ico_2_patch] = intersect(src(1).vertno, ico_2_labelno(i) + 1);       % left cortical patch
    c = 0;
    if isempty(ico_2_patch)                                                 
    [~,ico_2_patch] = intersect(src(2).vertno, ico_2_labelno(i) + 1);       % right cortical patch
    ico_2_patch = ico_2_patch + nuse;
    c = 1; 
    end
    ico(2).patch.no = ico_2_patch; n2 = [n2, ico_2_patch]; 
    if c == 0
        ico(2).labelno = src(1).vertno(ico(2).patch.no) - 1;
        ico(2).patch.neighbors = [];
        fwd_rr = fwd.src(1).rr; fwd_tris = fwd.src(1).tris;
        vertno = intersect(src(1).pinfo{ico(2).patch.no},fwd.src(1).vertno);
        rr = fwd.src(1).rr(vertno,:);
    else
        ico(2).labelno = src(2).vertno(ico(2).patch.no-nuse) - 1;
        fwd_rr = fwd.src(2).rr; fwd_tris = fwd.src(2).tris;
        vertno = intersect(src(2).pinfo{ico(2).patch.no-nuse},fwd.src(2).vertno);
        rr = fwd.src(2).rr(vertno,:);
    end
    ico(2).patch.neighbors = [];
    
    
    % Full Mesh
    if c == 0
        figure, set(gcf,'color','white');
        fig = gcf;
        trimesh(fwd_tris,fwd_rr(:,1),fwd_rr(:,2),fwd_rr(:,3),'FaceColor',[0.3 0.3 0.3],'EdgeColor','w');
        ax = fig.CurrentAxes;
        axis tight,  view(-90.5,2), camlight('left'); 
        hold all;  plot3(rr(:,1),rr(:,2),rr(:,3),'g.','MarkerSize',10); grid off;
    else
        figure, set(gcf,'color','white');
        fig = gcf;
        trimesh(fwd_tris,fwd_rr(:,1),fwd_rr(:,2),fwd_rr(:,3),'FaceColor',[0.3 0.3 0.3],'EdgeColor','w');
        ax = fig.CurrentAxes;
        axis tight%, view(ax, 90, 5, 2), camlight('right');
        hold all;  plot3(rr(:,1),rr(:,2),rr(:,3),'g.','MarkerSize',10); grid off;
    end
    
    
    pause; 
    clear ico_2_patch vertno rr ax fig rr fwd_rr
end

%% Save for Mask Creation
kk2 = [1 2 3 4];
save('icomasks_sim_patches.mat','ico','n2','kk2')