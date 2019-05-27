function angcfs_ov = cort_fscort_allpairs(sdiv, cpatch, savename)

%% Collect Eigenmodes for all Subcortical Volumes
Usfull = []; Ssfull = []; 
for j = 1:length(sdiv)
    if ~isempty(sdiv(j).U)
    Usfull = [Usfull sdiv(j).U];
    %Ssfull = [Ssfull sdiv(j).S]; %given different #eigenmodes difficult
    end
end

%% Angles Between Each Cortical Patch and All Subcortical Volumes
for i = 1:length(cpatch)
    disp(i);
    Uc = cpatch(i).U;
    Sc = cpatch(i).S;
    chk_excl = isempty(Uc);
    if chk_excl == 0
        angcfs_ov{i} = subspacea(Uc, Usfull)*180/pi;
    else
        angcc_ov{i,j} = [];
    end
    clear Uc Sc
end
save(savename,'angcfs_ov','-append');
clear Usfull %Ssfull

end