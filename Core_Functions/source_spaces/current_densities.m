function current_densities(sregname, datapath)

%% Current Densities By Region
cort_cd = 0.25*10^(-9);                                                     %[Am/mm^2] - surface density per Attal (2013)
hip_vcd = 0.25*10^(-9);                                                     %[Am/mm^3] - extrapolation
am_cd   = 0.25*10^(-9);                                                     %[Am/mm^3]
%am_cd   = 1.0*10^(-9);                                                     %[Am/mm^3]
put_cd  = 0.25*10^(-9);                                                     %[Am/mm^3]
caud_cd = 0.25*10^(-9);                                                     %[Am/mm^3] - addition
th_cd   = 0.025*10^(-9);                                                    %[Am/mm^3]
bs_cd   = 0.25*10^(-9);                                                     %[Am/mm^3] - addition
hip_scd = 1.0*10^(-9);                                                      %[Am/mm^2] - surface density per Matti and Okada (2016)
%hip_scd = 0.4*10^(-9);                                                     %[Am/mm^2] - surface density per Attal (2013)
lgn_cd  = 0.25*10^(-9);                                                     %[Am/mm^3] - ideally only principal orient
mgn_cd  = 0.25*10^(-9);                                                     %[Am/mm^3] - addition
ic_cd   = 0.25*10^(-9);                                                     %[Am/mm^3] - addition

%% Put all Volume Current Densities Together
for i = 1:length(sregname)
    if strcmp(sregname{i},'lh') || strcmp(sregname{i},'rh')
        subc_cd(i) = hip_vcd;
    elseif strcmp(sregname{i},'la') || strcmp(sregname{i},'ra')
        subc_cd(i) = am_cd;
    elseif strcmp(sregname{i},'lp') || strcmp(sregname{i},'rp')
        subc_cd(i) = put_cd; 
    elseif strcmp(sregname{i},'lc') || strcmp(sregname{i},'rc')
        subc_cd(i) = caud_cd; 
    elseif strcmp(sregname{i},'lt') || strcmp(sregname{i},'rt')
        subc_cd(i) = th_cd;
    elseif strcmp(sregname{i},'bs') || strcmp(sregname{i},'bsred')
        subc_cd(i) = bs_cd; 
    elseif strcmp(sregname{i},'lhipsurf') || strcmp(sregname{i},'rhipsurf')
        subc_cd(i) = hip_scd; 
    elseif strcmp(sregname{i},'llgn') || strcmp(sregname{i},'rlgn')
        subc_cd(i) = lgn_cd;
    elseif strcmp(sregname{i},'lmgn') || strcmp(sregname{i},'rmgn')
        subc_cd(i) = mgn_cd;
    elseif strcmp(sregname{i},'lic') || strcmp(sregname{i},'ric')
        subc_cd(i) = ic_cd;
    end
end
        
save(strcat(datapath,'/fwd/curr_den.mat'));

end
