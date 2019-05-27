function [U_subs_inds_reduced, flag_red] = reduce_combos(U_subs_inds, t_mode_inds)

for i = 1:length(U_subs_inds)
    mode_subset = U_subs_inds{i};
    for j = 1:size(t_mode_inds,1)
        %only include if a mode of any patch the top mode of the patch is included in the subset
        %s(i,j) = sum(ismember(t_mode_inds(j,:), mode_subset).*ismember(t_mode_inds(1,:), mode_subset)) > 0;
        for k = 1:size(t_mode_inds, 2)
            if ismember(t_mode_inds(j,k), mode_subset)                      %patch k mode j in mode_subset    
                s(j,k) = ismember(t_mode_inds(1,k), mode_subset);           %only accept if patch k mode 1 also included 
            else
                s(j,k) = 1;                                                 %not applicable, so dont make zero
            end
        end
    end
    flag_red(i) = prod(prod(s));
    %flag_red(i) = sum(s(i,:)); 
end

U_subs_inds_reduced = U_subs_inds(logical(flag_red > 0));
length(U_subs_inds_reduced)

end
      