function Xc_plot = srcspace_activity_surf(prefix, L, xcfull, cortp, cfwd, ico_num, source, ico, patch_exc, lkmax, fittype, Xs_source)

%% Initialize Parameters
patchno = (1:1:L);                                                          % all patch #s
patchno = patchno(patch_exc==0); L = length(patchno);                       % include all patch #s not on medial surface
trimP = 1:1:L;                                                              % Reg. MNE solution is not sparse
csrc = source(ico_num).hemisph;                                             % cortical source info structure
nuse = length(cfwd.src(1).vertno);                                          % # vertices in left hemisph
puse = length(csrc(1).pinfo);                                               % # patches in left hemisph
cnn = cat(1,csrc(1).nn,csrc(2).nn);                                         % all cortical orientations
%xcfull = repmat(xcfull,20);                                                % Repeat over many time points
T = size(xcfull,2);                                                         % # Time points - if T = 1 will generate a w file
xcfull = xcfull(1:lkmax,:);                                                 % Plot only cortical piece

%% Modes to Dipoles: Projection and Resultants
for ii = 1:L
    % Project Onto Modes
    V{ii} = ico(ico_num).patch(patchno(ii)).V(:,1:cortp);                   % Right eigenvector for chosen spatial modes
    Xc_plot(ii).mode = xcfull(cortp*(trimP(ii)-1)+1:cortp*trimP(ii));       % Mode Amplitudes
    Xc_plot(ii).dip = V{ii}*Xc_plot(ii).mode;                               % Dipole Amplitudes
    %Xc_plot(ii).dip_sc = Xc_plot(ii).dip/Xs_source;                        % ratio = cortical resultant fit/simulated subcortical source current
    
    % Compute Resultant for Patch
    Xc_plot(ii).orient = cnn(ico(ico_num).patch(patchno(ii)).lk,:);         % for orientation in MEG coordframe - use cfwd.source_nn(ico.patch(patch_num).lk,:) 
    nn = Xc_plot(ii).orient;                                                % orientation vector
    Xc_plot(ii).dip_or(:,:,1) = squeeze(nn(:,1)).*Xc_plot(ii).dip;          % vector of dipole current for orientation 1
    Xc_plot(ii).dip_or(:,:,2) = squeeze(nn(:,2)).*Xc_plot(ii).dip;          % vector of dipole current for orientation 2
    Xc_plot(ii).dip_or(:,:,3) = squeeze(nn(:,3)).*Xc_plot(ii).dip;          % vector of dipole current for orientation 3
    [Xc_plot(ii).dip_res, ~] = vec_sum(Xc_plot(ii).dip_or);                 % net dipole activity (amplitude and orientation) across all dipoles
    Xc_plot(ii).dip_res_sc = Xc_plot(ii).dip_res/Xs_source;                 % ratio = cortical resultant fit/simulated subcortical source current
end   
    
if T == 1                                                                   % Plot W File   
%% Create Vertex #s
w(1).vertices = cfwd.src(1).vertno;                                         % w vertex #s
w(2).vertices = cfwd.src(2).vertno;                                         % w vertex #s
    
%% Create W Data Fields
w(1).data = zeros(cfwd.src(1).nuse, 1);                                     % w data fields
w(2).data = zeros(cfwd.src(2).nuse, 1);                                     % w data fields
for ii = 1:L
	if patchno(ii) <= puse
    	w(1).data(ico(ico_num).patch(patchno(ii)).lk) = Xc_plot(ii).dip_res_sc;         %previously was dip_sc - i.e. each dipole to its own
    else
    	w(2).data(ico(ico_num).patch(patchno(ii)).lk-nuse) = Xc_plot(ii).dip_res_sc;    %previously was dip_sc - i.e. each dipole to its own
	end
end
    
%% Write W Files
mne_write_w_file1(strcat(prefix,'_ico_',num2str(ico_num),'_',fittype,'fit-lh.w'), w(1));
mne_write_w_file1(strcat(prefix,'_ico_',num2str(ico_num),'_',fittype,'fit-rh.w'), w(2));
else                                                                        % Plot STC File   
%% Initialize STC Parameters
stc(1).tmin = 0;  stc(2).tmin = stc(1).tmin;                                % seconds
stc(1).tstep = 0.01; stc(2).tstep = stc(1).tstep;                           % seconds
    
%% Create STC Vertex #s
stc(1).vertices = cfwd.src(1).vertno;                                       % stc vertex #s
stc(2).vertices = cfwd.src(2).vertno;                                       % stc vertex #s
    
%% Create STC Data Fields
stc(1).data = zeros(cfwd.src(1).nuse, T);                                   % stc data fields
stc(2).data = zeros(cfwd.src(2).nuse, T);                                   % stc data fields
for ii = 1:length(patchno)
	if patchno(ii) <= puse
    	stc(1).data(ico(ico_num).patch(patchno(ii)).lk,:) = Xc_plot(ii).dip_sc;
    else
    	stc(2).data(ico(ico_num).patch(patchno(ii)).lk-nuse,:) = Xc_plot(ii).dip_sc;
	end
end

%% Write STC Files
mne_write_stc_file1(strcat(prefix,'_ico_',num2str(ico_num),'_',fittype,'fit-lh.stc'), stc(1));
mne_write_stc_file1(strcat(prefix,'_ico_',num2str(ico_num),'_',fittype,'fit-rh.stc'), stc(2));
end

end