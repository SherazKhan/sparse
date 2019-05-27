function srcspace_figure(prefix, datatype, meastype, ico_num, append)

clc; close all;
dpath = strcat('/autofs/eris/purdongp/users/pavitrak/sourceloc_data/',datatype,'/',prefix);
add_on_name = []; %add_on_name = '_bs-full';                                %to define whether to use full or reduced brainstem
load(strcat(dpath,'/fwd/', prefix, '_subc_ico',num2str(ico_num),'_',meastype,'_prewhite_volumes',add_on_name,'.mat'));
T = 2; newn = zeros(256,256,256,T);                                         %T=2 even though we don't need time for coloring

for k = 1:length(sdiv_fwd)
    all_regions{k} = sdiv_fwd(k).reg_name;                                  %all regions stored in sdiv_fwd file
end

sregnums{1} = [2:6 8:12 14:15];                                             %IC only, and also Hipsurfs
savename{1} = ['srcspaces_setto_ico',num2str(ico_num),'_icam_hsurf',append];	   %IC only, and also Hipsurfs
if isempty(add_on_name)
sregnums{2} = [2:5 7:13];                                                   %reduced BS, and also Hipsurfs
savename{2} = ['srcspaces_setto_ico',num2str(ico_num),'_bsred_am_hsurf',append];   %reduced BS, and also Hipsurfs
else
sregnums{2} = [2:5 7:13];                                                   %full BS, and also Hipsurfs
savename{2} = ['srcspaces_setto_ico',num2str(ico_num),'_bsfull_am_hsurf',append];  %full BS, and also Hipsurfs
end

for i = 2:2
 sreg_case = sregnums{i};                                                   %the regions of interest for this case
 for sreg = 1:length(sreg_case)
  regname = sreg_fwd(sreg_case(sreg)).reg_name;                             %region of interest
  ind_sfwd = strncmp(regname, all_regions, 4);                              %indicate regions matching this one
  ssubdiv = find(ind_sfwd);                                                 %subdivisions of region i
  for j = 1:length(ssubdiv)
    maskname = strcat(regname,'_submask',num2str(j),'.mgz');
    A = MRIread(strcat(dpath,'/mri/svolume_decomps/',append,'/', maskname)); 
    M = A.vol == 1;                                                         %identify region where mask is high
    f = find(M==1);                                                         %find indices
    [x, y, z] = ind2sub(size(M),f);                                         %find x y z coordinates in voxel space
    
    for mm = 1:length(x)
        for t = 1:T
            a(t) = j*sreg+1;                                                %coloring should be different for each subdiv
            newn(x(mm),y(mm),z(mm),t) = a(t);                               %new mask
        end
    end
    clear M f x y z ssubdiv ind_sfwd
  end
 end
 B = A; B.vol = newn;                                                       %mask to write
 MRIwrite(B,strcat(dpath,'/mri/', savename{i},'.mgz'))                      %write mask in mgz format
end

end