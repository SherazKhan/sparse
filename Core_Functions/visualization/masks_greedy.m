function [mask, mask_offset, mask_mean] = masks_greedy(lk, clkmax, trimP, strimP, ctrimP, sdiv_fwd_plot, Xls_plot, ...
         mult, aseg, hlabel_fname, clabel_fname)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Derive coordinates of dipoles to create static masks of selected regions
% Written by Pavitra Krishnaswamy (pavitrak@nmr.mgh.harvard.edu)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Initialize Variables
srclocs = cell(1, length(trimP));                                           % initialize dipole coordinates in chosen regions
mask = cell(1, length(trimP));                                              % static mask for a chosen region
mask_offset = zeros(1,length(trimP)); mask_mean = mask_offset;              % offset and mean of time courses in a chosen region

for i = 1:length(trimP)
  s = strimP(i);                                                            % index wrt subcortical columns in joint srcspace
  c = ctrimP(i);                                                            % index wrt cortical columns in joint srscpace
 
  %% Derive Mask Offsets and Means
  mask_offset(i) = min(min(Xls_plot(i).dip_res));                           % no negative #s in masks - so need to shift
  mask_mean(i) = mean(mean(Xls_plot(i).dip_res));                           % mean value for scaling mask colors
  
  if s > 0                                                  
    %% Derive Static Masks for Subcortical Regions
    sdiv = sdiv_fwd_plot(s-sum(lk<=clkmax));                                % sdiv_fwd relevant to trimP(i)
    if isempty(strfind(sdiv.reg_name,'hipsurf'))                            % subcortical volumes
        srclocs{i} = sdiv.srclocs(1:3:end,1:3);                             % MRI coordinates (sRAS)
        sRAS = srclocs{i}';                                                 % put sRAS coordinates in correct format
        mricoords = mult(1:3,1:3)*sRAS;                                     % MRIcoords to voxels - without additive factor
        yyn = mricoords(1,:) + mult(1,4);                                   % MRIcoords to voxels - correct for coordinate frame
        xxn = mricoords(2,:) + mult(2,4);                                   % MRIcoords to voxels - correct for coordinate frame
        zzn = mricoords(3,:) + mult(3,4);                                   % MRIcoords to voxels - correct for coordinate frame
        newn = zeros(size(aseg.vol));                                       % initialize mask
        for mm = 1:length(xxn)                            
            newn(round(xxn(mm)), round(yyn(mm)), round(zzn(mm))) = 1;       % mask by index of dipole 
        end
        mask{i} = newn;                                                     % store subcortical vol. mask for future reference
    else                                                                    % hippocampal surfaces
        for hhpp_len = 1:length(hlabel_fname(i,:))
        maskname = strcat(pwd,'/labels/',hlabel_fname{i,hhpp_len}(1:strfind(hlabel_fname{i,hhpp_len},'.')-1),'.mgz');
        A = MRIread(maskname); 
        if hhpp_len == 1
            mask{i} = zeros(size(A.vol)); 
        end
        mask{i} = mask{i} + A.vol;                                          % store hippocampal mask for future reference
        end
    end
  else                                                                      % cortical case
    %% Derive Static Masks for Cortical Regions
    maskname = strcat(pwd,'/labels/',clabel_fname{i}(1:strfind(clabel_fname{i},'.')-1),'.mgz');
    A = MRIread(maskname); mask{i} = A.vol;                                 % store cortical mask for future reference
  end
  clear A maskname newn yyn xxn zzn sRAS mricoords;                         % clear all coordinates for next
end

end