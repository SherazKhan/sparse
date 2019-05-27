function [mask, mask_offset, mask_mean, vox_coord] =  masks_mne(lk, clkmax, trimP, strimP, ctrimP, sdiv_fwd_plot, ...
         Xls_plot, hsrc, hfwd, csrc, cfwd, mult, aseg)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Derive coordinates of dipoles to create static masks of selected regions
% Coordinate Frame Notes 
% src structures in mne fif files - coord_frame 5 = MRI coordinates = sRAS coordinates
% fwd structures in mne fif files - coord_frame 4 = MEG head coordinates
% references - pages 267 and 94 of MNE Manual 2.7.3
% Written by Pavitra Krishnaswamy (pavitrak@nmr.mgh.harvard.edu)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Initialize Variables
sRAS_coord = cell(1, length(trimP));                                        % initialize sRAS dipole coordinates 
vox_coord = cell(1, length(trimP));                                         % initialize voxel coordinates
mask = cell(1, length(trimP));                                              % static mask for a chosen region
mask_offset = zeros(1,length(trimP)); mask_mean = mask_offset;              % offset and mean of time courses in a chosen region
newn = zeros(size(aseg.vol));                                               % initialize mask

for i = 1:length(trimP)
  s = strimP(i);                                                            % index wrt subcortical columns in joint srcspace
  c = ctrimP(i);                                                            % index wrt cortical columns in joint srscpace
  
  %% Derive Mask Offsets and Means
  mask_offset(i) = min(min(Xls_plot(i).dip_res));                           % no negative #s in masks - so need to shift
  mask_mean(i) = mean(mean(Xls_plot(i).dip_res));                           % mean value for scaling mask colors
  
  if s > 0                                                  
    %% Derive sRAS Coordinates for Subcortical Regions 
    sdiv = sdiv_fwd_plot(s-sum(lk<=clkmax));                                % sdiv_fwd relevant to trimP(i)
    if isempty(strfind(sdiv.reg_name,'hipsurf'))                            % subcortical volumes
        %sRAS Coordinates
        sRAS_coord{i} = sdiv.srclocs(:,1:3)';                               % MRI coordinates (sRAS) put in correct format                                 
    else                                                                    % hippocampal surfaces 
        if isempty(strfind(sdiv.reg_name,'rh'));                            % left hippocampus
            if isfield(sdiv, 'patchno')
            [~,~,hlk] = intersect([hsrc(1).pinfo{sdiv.patchno}],hfwd.src(1).vertno,'legacy'); 
            else
            [~,~,hlk] = intersect(hsrc(1).pinfo{sdiv.index},hfwd.src(1).vertno,'legacy');
            end
            sRAS_coord{i} = hsrc(1).rr(hlk,:)'*1000;                        % MRI coordinates (sRAS for hippocampal patch)
        else                                                                % right hippocampus
            if isfield(sdiv, 'patchno')
            [~,~,hrk] = intersect([hsrc(2).pinfo{sdiv.patchno}],hfwd.src(2).vertno,'legacy'); 
            else
            [~,~,hrk] = intersect(hsrc(2).pinfo{sdiv.index},hfwd.src(2).vertno,'legacy');
            end
            sRAS_coord{i} = hsrc(2).rr(hrk,:)'*1000;                        % MRI coordinates (sRAS for hippocampal patch)
        end
    end
  else
     %% Derive sRAS Coordinates for Cortical Regions 
     patch_num = lk(c);                                                     % index of chosen patch
     nuse = length(csrc(1).pinfo);
     if patch_num <= nuse                                                   % left cortical patch
         [~,~,kk] = intersect(csrc(1).pinfo{patch_num},cfwd.src(1).vertno,'legacy');
         sRAS_coord{i} = csrc(1).rr(kk,:)'*1000;                            % MRI coordinates (sRAS for cortical patch)
     else                                                                   % right cortical patch
         [~,~,kk2] = intersect(csrc(2).pinfo{patch_num-nuse},cfwd.src(2).vertno,'legacy');
         sRAS_coord{i} = csrc(2).rr(kk2,:)'*1000;                           % MRI coordinates (sRAS for cortical patch)
     end
  end
  
  %% Voxel Coordinates
  mricoords = mult(1:3,1:3)*sRAS_coord{i};                                  % MRIcoords to voxels - without additive factor
  yyn = mricoords(1,:) + mult(1,4);                                         % MRIcoords to voxels - correct for coordinate frame
  xxn = mricoords(2,:) + mult(2,4);                                         % MRIcoords to voxels - correct for coordinate frame
  zzn = mricoords(3,:) + mult(3,4);                                         % MRIcoords to voxels - correct for coordinate frame
  vox_coord{i} = [xxn' yyn' zzn']';                                         % all voxel coordinates
  
  %% Static Mask
  for mm = 1:length(xxn)                            
      newn(round(xxn(mm)), round(yyn(mm)), round(zzn(mm))) = 1;             % mask by index of dipole 
  end
  clear yyn xxn zzn mricoords;                                              % clear all coordinates for next 
end
mask = newn;                                                                % store vol. mask for future reference

end