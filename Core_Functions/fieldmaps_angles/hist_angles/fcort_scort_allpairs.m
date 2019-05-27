function angfcs_ov = fcort_scort_allpairs(cpatch,sdiv, savename)

%% Collect Eigenmodes for all Cortical Patches
Ucfull = []; Scfull = []; 
for j = 1:length(cpatch)
    if ~isempty(cpatch(j).U)
	Ucfull = [Ucfull cpatch(j).U];
    %Scfull = [Scfull cpatch(j).S];
    end
end

%% Angles Between Each Subcortical Subdivision and All Cortical Patches     
for i = 1:length(sdiv)
	disp(i);
    Us = sdiv(i).U;
    Ss = sdiv(i).S;
    chk_excl = isempty(Us);
    if chk_excl == 0
        angfcs_ov{i} = subspacea(Us, Ucfull)*180/pi;
    else
        angfcs_ov{i} = [];
    end
    clear Us Ss
end
save(savename,'angfcs_ov','-append');
clear Ucfull %Scfull

end