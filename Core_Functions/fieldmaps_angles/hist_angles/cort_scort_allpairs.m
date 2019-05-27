function angcs_ov = cort_scort_allpairs(sdiv, cpatch, savename)

%% Angles Between Each Cortical Patch and Each Subcortical Subdivision
for i = 1:length(sdiv)
    for j = 1:length(cpatch)
        disp(i); 
        Us = sdiv(i).U;
        Ss = sdiv(i).S;
        disp(j);
        Uc = cpatch(j).U;
        Sc = cpatch(j).S;
        chk_excl = isempty(Uc) | isempty(Us);
        if chk_excl == 0
            angcs_ov{i,j} = subspacea(Us, Uc)*180/pi;
        else
            angcs_ov{i,j} = [];
        end
        clear Uc Sc Us Ss 
    end
end
save(savename,'angcs_ov','-append');

end