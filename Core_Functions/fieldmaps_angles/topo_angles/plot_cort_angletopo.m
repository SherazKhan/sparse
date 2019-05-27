function plot_cort_angletopo(source, ico, subdiv, patchno, fwd, prefix, w, cangle, flag)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% TOPOGRAPHIC MAP OF ANGLES
% Written by Pavitra Krishnaswamy and Gabriel Obregon-Henao
% pavitrak@nmr.mgh.harvard.edu, obregon@nmr.mgh.harvard.edu
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Compute Number of Patches in Left Hemisphere
nuse = length(fwd.src(1).vertno);
puse = length(source(subdiv).hemisph(1).pinfo);

%% Create W Vertex #s
w(1).vertices = fwd.src(1).vertno;
w(2).vertices = fwd.src(2).vertno;

%% Create W Data Fields
w(1).data = zeros(fwd.src(1).nuse, 1);
w(2).data = zeros(fwd.src(2).nuse, 1);
for ii = 1:length(patchno)
    if patchno(ii) <= puse
        w(1).data(ico(subdiv).patch(patchno(ii)).lk) = cangle(patchno(ii));    
    else
        w(2).data(ico(subdiv).patch(patchno(ii)).lk-nuse) = cangle(patchno(ii));
    end
end

%% Write W Files
mne_write_w_file1(strcat(prefix,'_ico_',num2str(subdiv),'_',flag,'-lh.w'), w(1));
mne_write_w_file1(strcat(prefix,'_ico_',num2str(subdiv),'_',flag,'-rh.w'), w(2));

end