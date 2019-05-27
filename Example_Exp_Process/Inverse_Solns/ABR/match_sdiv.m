clear; clc; close all; 
load('/autofs/eris/purdongp/users/pavitrak/sourceloc_data/AEP/cat_004/fwd/cat_004_subc_meg_eeg_plr-prewhite_SVDs_HFABR.mat');
for i = 1:length(sdiv_fwd)
    all_regions{i} = sdiv_fwd(i).reg_name;
end

ind_lic = find(strcmp(all_regions,'ic'),1);
ind_ric = find(strcmp(all_regions,'ic'),2); ind_ric = ind_ric(2); 
lic_locs = sdiv_fwd(ind_lic).srclocs(:,1:3);
lic_centroid = mean(lic_locs,1);
ric_locs = sdiv_fwd(ind_ric).srclocs(:,1:3);
ric_centroid = mean(ric_locs,1);

ind_bsred = find(strcmp(all_regions,'bsred'));
for i = 1:length(ind_bsred)
    curr_struct = sdiv_fwd(ind_bsred(i));
    bsdiv_locs = curr_struct.srclocs(:,1:3);
    bsdiv_centroid = mean(bsdiv_locs,1); 
    nl(i) = norm(bsdiv_centroid - lic_centroid);
    nr(i) = norm(bsdiv_centroid - ric_centroid);
end

[~,Il] = sort(nl);  %min will be Il(1)
[~,Ir] = sort(nr);  %min will be Il(2)

