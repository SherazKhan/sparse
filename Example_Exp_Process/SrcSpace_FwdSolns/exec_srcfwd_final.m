%% Initialization and Parameters
clear; clc; close all;
prefix = 'cat_004';                                                         %Subject name
meastype = 'meg_eeg';                                                       %MEG or EEG forward solution
datatype = 'AEP';                                                           %Paradigm of data: Illustrations, SEP, AEP
append = '_HFABR';                                                            %HFABR or MLR: indicator for diary name, meas_fname
ico_num = 3; hfwd_type = 'oct-3';                                           %cortical ico_num's that subcort srcspace division equalized to
%ico_num = 3; hfwd_type = 'ico-1';                                          %cortical ico_num's that subcort srcspace division equalized to
%ico_num = 2; hfwd_type = 'oct-1';                                          %cortical ico_num's that subcort srcspace division equalized to
est_ico_num = 3;                                                            %ico division - this is where estimation becomes 2 stage
whitener_fname = [];                                                        %non empty if whitener is stored separately than measurement file
add_allpaths;

if strcmp(append,'_HFABR')
    %% Cortical and Subcortical Source Space
    % Source and Fwd Files Creation
    auto_s_srcspace(prefix, datatype);                                      %Extract SrcSpace and Run Terminal Commands to Create SrcFwds
    read_save_sfwd(prefix, meastype, datatype, ico_num, append, hfwd_type);%Divide Deep Srcs
    if isempty(whitener_fname)
        svd_for_fwds(prefix, meastype, datatype, est_ico_num, append);      %use whitener and derive SVDs for cortical and subcortical fwds
    else
        svd_for_fwds(prefix, meastype, datatype, est_ico_num, append, whitener_fname); %use whitener and derive SVDs for cortical and subcortical fwds
    end
    
    % Visualizations and Summaries
    check_subdivs(prefix, datatype, meastype, ico_num, hfwd_type, append);  %MRI views for each subcortical subdivision and summary across deep sources
    srcspace_figure(prefix, datatype, meastype, ico_num, append);           %MRI views for all subdivided deep sources
else
    %% Cortical Source Space
    diary(strcat(prefix, '_cort_', meastype, '_srcspace_matlog',append));
    if isempty(whitener_fname)
        svd_for_fwds(prefix, meastype, datatype, est_ico_num, append);      %use whitener and derive SVDs for cortical and subcortical fwds
    else
        svd_for_fwds(prefix, meastype, datatype, est_ico_num, append, whitener_fname); %use whitener and derive SVDs for cortical and subcortical fwds
    end
    diary off;
end