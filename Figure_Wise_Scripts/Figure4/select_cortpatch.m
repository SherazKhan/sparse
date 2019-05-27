function lk = select_cortpatch(choose_lk, nc)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Select Cortical Patches
% Written Pavitra Krishnaswamy (pavitrak@nmr.mgh.harvard.edu)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

switch choose_lk
  
    case 'all'
    lk = 1:1:nc;                                                            %full iconum source space
  
    case 'sparse'
    path_patchest = '/autofs/eris/purdongp/users/pavitrak/sourceloc_scripts/figure_making/Figure3/som_erp_circ_results/';
    grp_fname = strcat(path_patchest,'SM04_forgrpangles_som_erp_meg_all_1g25m.mat');
    load(grp_fname, 'scortgrp','cortgrp','circ_name');
    lk = cortgrp(1).patchno;                                                %patch numbers for a sparse anatomical circuit
  
    case 'rand'
    lk = datasample(1:nc,5,'Replace',false);                                %patch numbers for a random set of 5 patches
    %datasample uses randperm to generate its values, hence not different techniques
    %http://www.mathworks.com/help/stats/datasample.html?requestedDomain=www.mathworks.com
  
    case 'none'
    lk = [];                                                                %no cortical patches at all
end

end