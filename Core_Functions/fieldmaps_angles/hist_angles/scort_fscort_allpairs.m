function angfss_ov = scort_fscort_allpairs(sdiv, savename)

%% Collect Eigenmodes for all Subcortical Volumes
Usfull = []; Ssfull = []; 
for j = 1:length(sdiv)
    if ~isempty(sdiv(j).U)
    Usfull = [Usfull sdiv(j).U];
    %Ssfull = [Ssfull sdiv(j).S]; %given different #eigenmodes difficult
    end
end

%% Angles Between Each Subcortical Subdivision and All Cortical Patches     
for i = 1:length(sdiv)
	disp(i);
    Us = sdiv(i).U;
    Ss = sdiv(i).S;
    ind = ones(1,length(sdiv));
    ind(i) = 0;
    chk_excl = isempty(Us);
    if chk_excl == 0
        angfss_ov{i} = subspacea(Us, Usfull(:,logical(ind)))*180/pi;
    else
        angfss_ov{i} = [];
    end
    clear Us Ss
end
save(savename,'angfss_ov','-append');
clear Usfull %Ssfull

end