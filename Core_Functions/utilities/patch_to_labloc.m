function [lblno, vertno, srclocs, lr] = patch_to_labloc(src, fwd, patch_num)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Given a patch number corresponding to srcspace (src) and fwd structure
% (fwd), identify label numbers, vertex numbers and sRAS source locations
% Written by Gabriel Obregon and Pavitra Krishnaswamy
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

puse = length(src(1).pinfo);
if patch_num <= puse
    lr = 'l';                                                                   %left
    lblno = src(1).vertno(patch_num)-1;                                         %label number
    vertno = intersect(src(1).pinfo{patch_num},fwd.src(1).vertno,'legacy');     %vertex number
    srclocs =  fwd.src(1).rr(vertno,:);                                         %MRI coordinates
else
    lr = 'r';                                                                   %right
    lblno = src(2).vertno(patch_num-puse)-1;                                    %label number
    vertno = intersect(src(2).pinfo{patch_num-puse},fwd.src(2).vertno,'legacy');%vertex number
    srclocs =  fwd.src(2).rr(vertno,:);                                         %MRI coordinates
end

end