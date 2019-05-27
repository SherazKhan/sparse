function patchno = lab_to_patch(src, lblno, lr)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Given a label number corresponding to srcspace (src) identify patch no. 
% src(1).vertno(patch_num) should be lblno + 1
% Written by Pavitra Krishnaswamy
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

puse = length(src(1).pinfo); 
for i = 1:length(lblno)
  if lr(i) == 1
    [~, ind, ~] = intersect(src(1).vertno, lblno(i)+1,'legacy');
    patchno(i) = ind;
  else
    [~, ind, ~] = intersect(src(2).vertno, lblno(i)+1,'legacy');
    patchno(i) = ind + puse;
  end
end

end