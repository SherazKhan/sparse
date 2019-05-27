function [c_currstr, s_currstr] = currstr_calc(datapath, prefix, meastype, chtype, gradmag, ico_num, sregname, sregsdiv)

load(strcat(datapath, '/fwd/curr_den.mat'),'cort_cd', 'th_cd');             % Read in cortical and VPL current densities

load(strcat(datapath, '/fwd/', prefix,'_cort_meg_',chtype,'_',gradmag,'_SVDs.mat'),'ico'); % Read in Cortical forward solutions
c_currstr = cort_cd*mean([ico(ico_num).patch(:).A]);                        % mean cortical current strength across patches 

load(strcat(datapath, '/fwd/', '2014-15_Vectorview_System/whiten_',gradmag,...                                        
    '/',prefix,'_subc_',meastype,'_',chtype,'_', gradmag,'_SVDs.mat'), 'sdiv_fwd'); % 2014 Thesis Version for Subcortical SVDs
for i = 1:length(sdiv_fwd) 
    all_regions{i} = sdiv_fwd (i).reg_name;                                 % All Region Names
end
ind = find(strcmp(sregname, all_regions));                                  % Identify all region names
sdiv = ind(sregsdiv); clear ind;                                            % Identify subdivision index amongst sdiv_fwd divisions
s_currstr = th_cd*sdiv_fwd(sdiv).volume;                                    % VPL Current Strength

end