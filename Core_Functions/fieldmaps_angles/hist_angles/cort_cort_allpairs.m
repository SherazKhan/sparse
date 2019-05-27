function angcc_ov = cort_cort_allpairs(cpatch, savename)

%% Angles Between Each Cortical Patch and Every Other Cortical Patch
for i = 1:length(cpatch)
    for j = 1:length(cpatch) 
        disp(i); 
        Uc1 = cpatch(i).U;
        Sc1 = cpatch(i).S;
        disp(j);
        Uc2 = cpatch(j).U;
        Sc2 = cpatch(j).S;
        chk_excl = isempty(Uc1) | isempty(Uc2);
        if chk_excl == 0
            angcc_ov{i,j} = subspacea(Uc1, Uc2)*180/pi;
        else
            angcc_ov{i,j} = [];
        end
        clear Uc1 Sc1 Uc2 Sc2
    end
end
save(savename,'angcc_ov','-append');

end