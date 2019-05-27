function [div_subs, dflag] = all_combs(divs_struct)

l = length(divs_struct);                                                    %number of divisions
reg_inds = 1:1:l;                                                           %indices of divisions
cnt = 1;

%% Combinations For Comparisons with Other Div_Struct
for ll = 1:l
    combs_possible{ll} = nchoosek(reg_inds,ll);
    for i = 1:size(combs_possible{ll},1)
  	  divs_sel{cnt} = combs_possible{ll}(i,:);                          %present subset of divisions
      div_subs{cnt} = divs_struct(divs_sel{cnt});                           %structure for present subset of divisions
      cnt = cnt + 1;
    end   
end
 
%% Flag for Combinations For Comparisons with Itself
for i = 1:cnt-1
    for j = 1:cnt-1
    %flag if this combination can be used
    dflag(i,j) = sum(ismember(divs_sel{i}, divs_sel{j})) > 0; 
    end
end

end