function [all_regions, sfwd, sp] = collate_sfwd_storage(sfwd, raw, varargin)
  num = length(sfwd); all_regions = cell(1, num); 
  regvec = {sfwd.reg_name}; 
  
  if nargin == 2
    if raw == 0
      sp = [sfwd.numwhite_modes];                               %#white modes for 95% NMRA
    else
      sp = [sfwd.numraw_modes];                                 %#raw modes for 95% NMRA
    end
  else
      sp = repmat(varargin{1},1, num);                          %# modes defined by input 
  end
  init = zeros(max(sp),max(sp));                                %initialize singular vals matrix
  
  for i = 1:num
      if ~isempty(regvec{i})                                    % if fwd for region is in sfwd
      all_regions{i} = sfwd(i).reg_name;                        %name of region of interest
      if raw == 0
          sfwd(i).whiteUs_p = sfwd(i).whiteUs(:,1:sp(i));       %white singular vectors
          sfwd(i).whiteS_p = init;                              %initialize singular vals matrix
          sfwd(i).whiteS_p(1:sp(i),1:sp(i)) = sfwd(i).whiteS(1:sp(i),1:sp(i));   %white singular values
          sfwd(i).whiteGstheta_p = sfwd(i).whiteGstheta(:,1:sp(i));  %reduced rank approx to white cort fwd
      else
          sfwd(i).rawUs_p = sfwd(i).rawUs(:,1:sp(i));           %raw singular vectors
          sfwd(i).rawS_p = init;                                %initialize singular vals matrix
          sfwd(i).rawS_p(1:sp(i),1:sp(i)) = sfwd(i).rawS(1:sp(i),1:sp(i));   %raw singular values
          sfwd(i).rawGstheta_p = sfwd(i).rawGstheta(:,1:sp(i)); %reduced rank approx to raw cort fwd
      end
      end
  end
end