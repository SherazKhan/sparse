function [angs_allcombos_sc, angs_allcombos_cc] = compute_angs_allcombos(scortgrp, cortgrp, type_flag)

%% Read Subcortical and Cortical Components of Circuit and Regions to Analyze
sdiv = scortgrp.sdiv;                                                       %subcortical divisions in circuit
scort_p = sdiv.sp;                                                          %number of eigenmodes in a division
cpatch = cortgrp.cpatch;                                                    %cortical patches in circuit        
cort_p = cpatch(1).cp;                                                      %number of eigenmodes in a patch
%if ~strcmp(type_flag, 'all') && ~strcmp(type_flag,'top') && ~strcmp(type_flag,'rand') %each patch or division to be represented by all its modes
%    [sdiv_subs, sflag] = all_combs(sdiv);                                  %all subsets of subcortical divisions
%    [cdiv_subs, cflag] = all_combs(cpatch);                                %all subsets of subcortical divisions
%end

%% Angles for All Combinations of Modes in Source Space: Scort vs. Cort
disp('subcort vs cort'); 
angs_allcombos_sc = reg1_reg2_allcombos(sdiv, cpatch, type_flag, scort_p, cort_p);%angles for all combinations of modes

%% Angles for All Combinations of Modes in Source Space: Cort vs. Cort Control
disp('cort vs cort'); 
angs_allcombos_cc = reg1_reg2_allcombos(cpatch, cpatch, type_flag, cort_p, cort_p); %angles for all combinations of modes

end