function sras_coord = reg_coords(lk_flag, lk, clkmax, sdiv_fwd, hsrc, hfwd, csrc, cfwd)

%[lk2, lkind, lk2ind] = unique(lk); lk_flag = lk_flag(lkind);                      %find unique values of lk and lk_flag by division
%if use the above, then one must use: sras_coord = cell(1, length(lk); and assign sras_coord{lk2ind} to the sras_coord we calculate below

for i = 1:length(lk)    
  if ~lk_flag(i)
    %% Derive sRAS Coordinates for Subcortical Regions 
    sdiv = sdiv_fwd(lk(i)-clkmax);                                          % sdiv_fwd relevant to trimP(i)
    if isempty(strfind(sdiv.reg_name,'hipsurf'))                            % subcortical volumes
        %sRAS Coordinates
        sras_coord{i} = sdiv.srclocs(1:3:end,1:3)';                         % MRI coordinates (sRAS) put in correct format                                 
    else                                                                    % hippocampal surfaces 
        if isempty(strfind(sdiv.reg_name,'rh'));                            % left hippocampus
            [~,~,hlk] = intersect([hsrc(1).pinfo{[sdiv.index]}],hfwd.src(1).vertno,'legacy'); 
            sras_coord{i} = hsrc(1).rr(hlk,:)'*1000;                        % MRI coordinates (sRAS for hippocampal patch)
        else                                                                % right hippocampus
            [~,~,hrk] = intersect([hsrc(2).pinfo{[sdiv.index]}],hfwd.src(2).vertno,'legacy'); 
            sras_coord{i} = hsrc(2).rr(hrk,:)'*1000;                        % MRI coordinates (sRAS for hippocampal patch)
        end
    end
  else
     %% Derive sRAS Coordinates for Cortical Regions 
     patch_num = lk(i);                                                     % index of chosen patch
     nuse = length(csrc(1).pinfo);
     if patch_num <= nuse                                                   % left cortical patch
         [~,~,kk] = intersect(csrc(1).pinfo{patch_num}, cfwd.src(1).vertno, 'legacy');
         sras_coord{i} = csrc(1).rr(kk,:)'*1000;                            % MRI coordinates (sRAS for cortical patch)
     else                                                                   % right cortical patch
         [~,~,kk2] = intersect(csrc(2).pinfo{patch_num-nuse}, cfwd.src(2).vertno, 'legacy');
         sras_coord{i} = csrc(2).rr(kk2,:)'*1000;                           % MRI coordinates (sRAS for cortical patch)
     end
  end
end
 
end