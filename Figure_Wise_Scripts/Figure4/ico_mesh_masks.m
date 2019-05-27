%% Written by Gabriel Obregon-Henao (obregon@nmr.mgh.harvard.edu)
close all; clear all; clc
subjid = 'SM04';
dir_path = strcat('/autofs/cluster/purdonlab/users/pavitrak/sourceloc_data/Illustrations/',subjid);
addpath(genpath(dir_path))
src2_fname = strcat(subjid,'-ico-2p-src.fif');         csrc2 = mne_read_source_spaces(src2_fname);
fwd_fname = strcat(subjid,'-meg-all-fixed-fwd.fif');   cfwd = mne_read_forward_solution(fwd_fname);
aseg = MRIread(strcat(dir_path,'/mri/aseg.mgz'));                           % all masks information
load(strcat('T3.mat'));                                                     % sRAS to vRAS transform
v2ras1 = aseg.vox2ras1; mult = inv(v2ras1)*T3;                              % for voxel coordinate conversion       
load icomasks_sim_patches;                                                  % file containing patch # information across ico's

for ico_num = 2:2%1:3
%% Switch by ICO
switch ico_num
    case 1
        patchno = n1; csrc = csrc1; sel = kk1; lab = 'ico1';                % index of chosen patch - n1(kk1)
    case 2
        patchno = n2; csrc = csrc2; sel = kk2; lab = 'ico2';                % index of chosen patch - n2(kk2)
    case 3
        patchno = n3; csrc = csrc3; sel = kk3; lab = 'ico3';                % index of chosen patch - n3(kk3)
end

for i = 1:length(patchno)
patch_num = patchno(i);                                                     % patch # of interest
newn = zeros(size(aseg.vol));                                               % initialize newn
%% Derive sRAS Coordinates for Cortical Regions
nuse = length(csrc(1).pinfo);
if patch_num <= nuse                                                        % left cortical patch
    [~,~,kkl] = intersect(csrc(1).pinfo{patch_num},cfwd.src(1).vertno);
    sRAS_coord{i} = csrc(1).rr(kkl,:)'*1000;                                % MRI coordinates (sRAS for cortical patch)
else                                                                        % right cortical patch
	[~,~,kkr] = intersect(csrc(2).pinfo{patch_num-nuse},cfwd.src(2).vertno);
    sRAS_coord{i} = csrc(2).rr(kkr,:)'*1000;                                % MRI coordinates (sRAS for cortical patch)
end
     
%% Voxel Coordinates
mricoords = mult(1:3,1:3)*sRAS_coord{i};                                    % MRIcoords to voxels - without additive factor
yyn = mricoords(1,:) + mult(1,4);                                           % MRIcoords to voxels - correct for coordinate frame
xxn = mricoords(2,:) + mult(2,4);                                           % MRIcoords to voxels - correct for coordinate frame
zzn = mricoords(3,:) + mult(3,4);                                           % MRIcoords to voxels - correct for coordinate frame
vox_coord{i} = [xxn' yyn' zzn']';                                           % all voxel coordinates
  
%% Static Mask
for mm = 1:length(xxn)                            
    newn(round(xxn(mm)), round(yyn(mm)), round(zzn(mm))) = 1;               % mask by index of dipole 
end
B = aseg; B.vol = newn;
if sel == i
    savename = strcat(lab,'masksel',num2str(patch_num),'.mgz');             % selected mask - color differently
else
    savename = strcat(lab,'mask',num2str(patch_num),'.mgz');                % random mask
end
MRIwrite(B,savename);                                                       % write MRI mask for labels of interest

clear mricoords yyn xxn zzn newn B
end
clear sRAS_coord vox_coord 
end