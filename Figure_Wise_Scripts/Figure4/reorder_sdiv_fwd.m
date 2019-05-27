%% Paths and Parameters
clear; clc; close all;
prefix = 'SM04'; meastype = 'meg'; datatype = 'Illustrations';
g = 1; m = 25; gradmag = strcat(num2str(g), 'g', num2str(m), 'm'); chtype = 'all'; %scaling ratios for gradiom vs. magnetom                                                            
dpath = '/autofs/eris/purdongp/users/pavitrak/'; add_allpaths;              %add all paths
fwd_path = strcat(datapath, '/fwd/');                                       %latest forward path

%% Load Sdiv_Fwd File 
ssvds_fname = strcat(fwd_path, prefix, '_subc_ico3_meg_', chtype, '_', gradmag, '_SVDs.mat');
load(ssvds_fname,'sdiv_fwd','sregname');                                    %subcortical fwds by subdivisions with regnames
for k = 1:length(sdiv_fwd)
    all_regions{k} = sdiv_fwd(k).reg_name;                                  %save all region names in detail
end
    
%% Rearrange Regions 
% currently: la lp lc lt ra rp rc rt bsred lhipsurf rhipsurf
% rearrange: lhipsurf la lp lc lt bsred rt rc rp ra rhipsurf
reordered_sregname = {'lhipsurf','la','lp','lc','lt','bsred','rt','rc','rp','ra','rhipsurf'};
cnt = 1; 
for i = 1:length(reordered_sregname)
    regname{i} = reordered_sregname{i};                                     %name of selected region
    ind_sfwd = strncmp(regname{i}, all_regions, 4);                         %find indices of sdiv_fwd that match regname{i}
    ind = find(ind_sfwd);                                                   %indices in sdiv_fwd corresponding to this region
    if strcmp(regname{i},'bsred')
        sdiv = sdiv_fwd(ind);
        for j = 1:length(sdiv)
            srcloc_bs{j} = [sdiv(j).srclocs(:,1)];                          %L-R axis is 1st srclocs coordinate
            lr_bs(j) = sum([srcloc_bs{j}(:)] < 1) > round(0.5*length(srcloc_bs{j})); %L negative, R positive, lr_bs is 1 if majority dipoles left
        end
        sdiv_bs = [sdiv(lr_bs) sdiv(~lr_bs)];                               %left first right next
    end
    cnt_st = cnt; cnt_en = cnt_st + length(ind) - 1;                        %indices in new sdiv_fwd
    if strcmp(regname{i},'bsred')
        reordered_sdiv_fwd(cnt_st:cnt_en) = sdiv_bs;                        %Store brainstem
    else
        reordered_sdiv_fwd(cnt_st:cnt_en) = sdiv_fwd(ind);                  %Store ith Region
    end
    cnt = cnt_en + 1;
end
clear sdiv_fwd sregname;

%% Save Reordered Sdiv_Fwd
sdiv_fwd = reordered_sdiv_fwd;
sregname = reordered_sregname;
savename = strcat(fwd_path, prefix, '_subc_ico3_meg_', chtype, '_', gradmag, '_SVDs_rearranged.mat');
save(savename,'sdiv_fwd','sregname');
for k = 1:length(sdiv_fwd)
    reordered_all_regions{k} = sdiv_fwd(k).reg_name;
end