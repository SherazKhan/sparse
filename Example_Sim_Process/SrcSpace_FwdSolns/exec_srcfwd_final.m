%% Initialization and Parameters
clear; clc; close all;
prefix = 'SM04';                                                            %Subject name
meastype = 'meg';                                                           %MEG or EEG forward solution
datatype = 'SEP';                                                           %Paradigm of data: Illustrations, SEP, AEP
overlap = '_overlap';    %overlap = '_overlap' or '_no-overlap';            %indicator for diary name, meas_fname (read_save_fwd, cort_surf_fwd)
add_allpaths;
ico_num = 3; hfwd_type = 'ico-1';                                           %cortical ico_num's that subcort srcspace division equalized to
%ico_num = 2; hfwd_type = 'oct-1';                                          %cortical ico_num's that subcort srcspace division equalized to
est_ico_num = 3;                                                            %ico division - this is where estimation becomes 2 stage
whitener_fname = [];                                                        %non empty if whitener is stored separately than measurement file

%% Source and Fwd Files Creation
nuc_srcspace(prefix, datatype);                                             %IC
auto_s_srcspace(prefix, datatype);                                          %Extract SrcSpace and Run Terminal Commands to Create SrcFwds
read_save_sfwd(prefix, meastype, datatype, ico_num, overlap, hfwd_type);    %Divide Deep Srcs - !!update meas_filename in lines 60-62 of read_save_sfwd
if isempty(whitener_fname)
    svd_for_fwds(prefix, meastype, datatype, est_ico_num, overlap);         %use whitener and derive SVDs for cortical and subcortical fwds
else
    svd_for_fwds(prefix, meastype, datatype, est_ico_num, overlap, whitener_fname); %use whitener and derive SVDs for cortical and subcortical fwds
end

%% Visualizations and Summaries
check_subdivs(prefix, datatype, meastype, ico_num, hfwd_type, overlap);     %MRI views for each subcortical subdivision and summary across deep sources
srcspace_figure(prefix, datatype, meastype, ico_num, overlap);              %MRI views for all subdivided deep sources