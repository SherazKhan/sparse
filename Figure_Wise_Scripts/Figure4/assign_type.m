function [resmat_name, cax, dist_name, savefigname, savestatsname] = assign_type(testcase, prefix, savepath)

switch testcase
    case 1                                                                  % All Cortical and All Subcortical
    choose_cort = 'all';  choose_scort = 'all';            
    resmat_name = strcat(savepath, prefix, '_', choose_cort, '_cort_', choose_scort, '_scort_resmat.mat'); 
    dist_name = resmat_name;
    cax = [0 0.035]; %cax = [0 0.007];%[-0.0025 0.01]; %  
    
    case 2                                                                  % Sparse Cortical and All Subcortical
    choose_cort = 'sparse';   choose_scort = 'all';    
    resmat_name = strcat(savepath, prefix, '_', choose_cort, '_cort_', choose_scort, '_scort_resmat.mat'); 
    dist_name = resmat_name;
    cax = [0 0.1];%cax = [0 0.07];%[-0.02 0.1]; %  
    
    case 3                                                                  % Greedy Sparse Cortical and Sparse Subcortical
    choose_cort = 'sparse';   choose_scort = 'sparse';    
    resmat_name = strcat(savepath, prefix, '_',choose_cort, '_cort_', choose_scort,'_scort_resmat.mat');
    dist_name = strcat(savepath, prefix, '_',choose_cort, '_cort_all_scort_resmat.mat');
	cax = [0 1];%[-0.05 1];%
    
    case 4                                                                  % Cortical only
    choose_cort = 'all'; choose_scort = 'none';            
    resmat_name = strcat(savepath, prefix, '_', choose_cort, '_cort_', choose_scort, '_scort_resmat.mat'); 
    dist_name = resmat_name;
    cax = [0 0.15]; 
    
    case 5                                                                  % Subcortical only
    choose_cort = 'none'; choose_scort = 'all';             
    resmat_name = strcat(savepath, prefix, '_', choose_cort, '_cort_', choose_scort, '_scort_resmat.mat'); 
    dist_name = resmat_name;
    cax = [0 0.1];
end
savefigname =   strcat(savepath, prefix, '_', choose_cort, '_cort_', choose_scort, '_scort_K');
savestatsname = strcat(savepath, prefix, '_', choose_cort, '_cort_', choose_scort, '_scort_summ');
    
end