function trimP_Inds = id_repmodes(reg_index, lk, svolume, cort_scort, corrind)
%function [trimP_Inds, Subdiv_Inds, Chosen_Inds, Repeat_Inds] = id_repmodes(reg_index, lk, svolume, cort_scort, corrind)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Written by Pavitra Krishnaswamy (pavitrak@nmr.mgh.harvard.edu)
% Reads reg_index, chosen lk, and the svolume -> to identify chosen modes from repeating subdivision
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[~,Repeat,ic] = unique(reg_index, 'rows');
countofx = hist(ic, 1:max(ic));
Ind_Repeat_Rows = (countofx ~= 1);
Repeat_Rows = Repeat(Ind_Repeat_Rows);
Repeat_Counts = countofx(Ind_Repeat_Rows);
[~, ind2] = find(corrind > 0);

if isempty(Repeat_Rows)
    trimP_Inds = [];
else
for i = 1:length(Repeat_Rows)
	ind = find(ic == ic(Repeat_Rows(i)));
    Chosen_Inds{i} = lk(ind);
    if strcmp(cort_scort,'s')
        [~, Repeat_Inds{i}, ~] = intersect(Chosen_Inds{i}, svolume{reg_index(ind(1),1), reg_index(ind(1),2)}, 'stable');%,'legacy'); 
        %changed from stable and legacy on 21 Dec 2015
        Subdiv_Inds{i} = Chosen_Inds{i}(Repeat_Inds{i});
    else
        Repeat_Inds{i} = Chosen_Inds{i};
        Subdiv_Inds{i} = Chosen_Inds{i};
    end
    ind = ind2(ind);
    trimP_Inds{i} = ind; clear ind
end
end

end

