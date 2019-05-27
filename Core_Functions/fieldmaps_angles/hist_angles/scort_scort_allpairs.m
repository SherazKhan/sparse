function angss_ov = scort_scort_allpairs(sdiv, savename)

%% Angles between All Pairs of Subcortical Divisions
for i = 1:length(sdiv)
    for j = 1:length(sdiv)
        disp(i); 
        Us1 = sdiv(i).U;
        Ss1 = sdiv(i).S;
        disp(j);
        Us2 = sdiv(j).U;
        Ss2 = sdiv(j).S;
        chk_excl = isempty(Us1) | isempty(Us2);
        if chk_excl == 0
            angss_ov{i,j} = subspacea(Us1, Us2)*180/pi;
        else
            angss_ov{i,j} = [];
        end
        clear Us1 Ss1 Us2 Ss2
    end
end
save(savename,'angss_ov','-append');

end