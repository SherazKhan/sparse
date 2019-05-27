function Xc_plot = srcspace_activity_mri(prefix, L, xcfull, cortp, cfwd, ico_num, source, ico, mri_path, patch_exc, lkmax, fittype, Xs_source)

%% Initialize Parameters
patchno = (1:1:L);                                                          % all patch #s
patchno = patchno(patch_exc==0); L = length(patchno);                       % include all patch #s not on medial surface
trimP = 1:1:L;                                                              % Reg. MNE solution is not sparse
csrc = source(ico_num).hemisph;                                             % cortical source info structure
nuse = length(cfwd.src(1).vertno);                                          % # vertices in left hemisph
puse = length(csrc(1).pinfo);                                               % # patches in left hemisph
cnn = cat(1,csrc(1).nn,csrc(2).nn);                                         % all cortical orientations
ico = ico(ico_num);                                                         % consider only ico-num
aseg = MRIread(strcat(mri_path,'aseg.mgz'));                                % aseg for mask creation
load(strcat(mri_path,'T3.mat'));                                            % sRAS to vRAS transform
v2ras1 = aseg.vox2ras1; mult = inv(v2ras1)*T3;                              % transform between sRAS and voxels
Xc_plot(1:L) = struct('mode',[],'orient',[],'dip',[],'dip_res',[]);         % Initialize Structure
xcfull = xcfull(1:lkmax,:);                                                 % Plot only cortical piece

%% Modes to Dipoles: Projection and Resultants
for ii = 1:L
    % Project Onto Modes
    patch_num = patchno(ii);
    V{ii} = ico.patch(patch_num).V(:,1:cortp);                              % Right eigenvector for chosen spatial modes
    Xc_plot(ii).mode = xcfull(cortp*(trimP(ii)-1)+1:cortp*trimP(ii));       % Mode Amplitudes
    Xc_plot(ii).dip = V{ii}*Xc_plot(ii).mode;                               % Dipole Amplitudes
    
    % Compute Resultant for Patch
    Xc_plot(ii).orient = cnn(ico.patch(patch_num).lk,:);                    % for orientation in MEG coordframe - use cfwd.source_nn(ico.patch(patch_num).lk,:) 
    nn = Xc_plot(ii).orient;                                                % orientation vector
    Xc_plot(ii).dip_or(:,:,1) = squeeze(nn(:,1)).*Xc_plot(ii).dip;          % vector of dipole current for orientation 1
    Xc_plot(ii).dip_or(:,:,2) = squeeze(nn(:,2)).*Xc_plot(ii).dip;          % vector of dipole current for orientation 2
    Xc_plot(ii).dip_or(:,:,3) = squeeze(nn(:,3)).*Xc_plot(ii).dip;          % vector of dipole current for orientation 3
    [Xc_plot(ii).dip_res, ~] = vec_sum(Xc_plot(ii).dip_or);                 % net dipole activity (amplitude and orientation) across all dipoles
    Xc_plot(ii).dip_res_sc = Xc_plot(ii).dip_res/Xs_source;                 % ratio = cortical resultant fit/simulated subcortical source current
end

for ii = 1:L
%% Derive sRAS Coordinates for Cortical Regions 
patch_num = patchno(ii);
nuse = length(csrc(1).pinfo);
if patch_num <= nuse                                                        % left cortical patch
	[~,~,kk] = intersect(csrc(1).pinfo{patch_num}, cfwd.src(1).vertno,'legacy');
    sRAS_coord{ii} = csrc(1).rr(kk,:)'*1000;                                % MRI coordinates (sRAS for cortical patch)
else                                                                        % right cortical patch
	[~,~,kk2] = intersect(csrc(2).pinfo{patch_num-nuse}, cfwd.src(2).vertno,'legacy'); 
    sRAS_coord{ii} = csrc(2).rr(kk2,:)'*1000;                               % MRI coordinates (sRAS for cortical patch)
end
   
%% Voxel Coordinates
mricoords = mult(1:3,1:3)*sRAS_coord{ii};                                   % MRIcoords to voxels - without additive factor
yyn = mricoords(1,:) + mult(1,4);                                           % MRIcoords to voxels - correct for coordinate frame
xxn = mricoords(2,:) + mult(2,4);                                           % MRIcoords to voxels - correct for coordinate frame
zzn = mricoords(3,:) + mult(3,4);                                           % MRIcoords to voxels - correct for coordinate frame
vox_coord{ii} = [xxn' yyn' zzn']';                                          % all voxel coordinates

%% Static Mask
for mm = 1:length(xxn)
	mrimask(round(xxn(mm)), round(yyn(mm)), round(zzn(mm))) = Xc_plot(ii).dip_res_sc; % mask by index of dipole 
end
clear yyn xxn zzn mricoords;                                                % clear all coordinates for next patch
end

%% Save Mask Movie
B = aseg; B.vol = mrimask;                                                  % Initialize Mask
savename = strcat(prefix, '_ico_',num2str(ico_num),'_',fittype,'fit_mask.mgz'); % Name of MRI Mask File for Cortical Fit
MRIwrite(B,savename);                                                       % Write MRI mask for Displaying Cortical Fit
disp('Wrote Mask Movie to Display MNE Estimates');                          % Completed

end