function angfcc_ov = fcort_cort_allpairs(cpatch, savename)

%% Collect Eigenmodes for all Cortical Patches
Ucfull = []; Scfull = []; 
for j = 1:length(cpatch)
    if ~isempty(cpatch(j).U)
	Ucfull = [Ucfull cpatch(j).U];
    %Scfull = [Scfull cpatch(j).S];
    end
end

%% Angles Between Each Subcortical Subdivision and All Cortical Patches     
for i = 1:length(cpatch)
	disp(i);
    Uc = cpatch(i).U;
    Sc = cpatch(i).S;
    ind = ones(1,length(cpatch));
    ind(i) = 0;
    chk_excl = isempty(Uc);
    if chk_excl == 0
        angfcc_ov{i} = subspacea(Uc, Ucfull(:,logical(ind)))*180/pi;
    else
        angfcc_ov{i} = [];
    end
    clear Uc Sc
end
save(savename,'angfcc_ov','-append');
clear Ucfull %Scfull

end