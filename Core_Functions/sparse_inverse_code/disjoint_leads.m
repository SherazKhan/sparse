function [Gn,Gtheta,KK,currstr] = disjoint_leads(ico,p,pt_s,patchno,source,subdiv,subdiv_sol, curr_den_fname)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% FORM REDUCED-ORDER LEAD FIELD OF CANDIDATE SUPPORT OVER SUBDIV 
% Written by Gabriel Obregon-Henao (obregon@nmr.mgh.harvard.edu)
% Updated by Pavitra Krishnaswamy
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Initialize
pk = [];
load(curr_den_fname, 'cort_cd');                                            %[Am/mm^2] - cortical current density

%% Form Candidate Support
if isequal(subdiv,1) 
    % Load candidate support in ico-1p
    pk = 1:length(ico(subdiv).patch);      
else
    % Compute number of patches in left hemisphere of subdiv-1
    puse = length(source(subdiv-1).hemisph(1).pinfo);
    
    % Find all patches in subdiv contained in estimate over subdiv-1
    for ii = 1:length(patchno)
        if patchno(ii) <= puse
            [~,~,lk] = intersect(source(subdiv-1).hemisph(1).pinfo{patchno(ii)},source(subdiv).hemisph(1).vertno,'legacy');
        else
            [~,~,lk] = intersect(source(subdiv-1).hemisph(2).pinfo{patchno(ii)-puse},source(subdiv).hemisph(2).vertno,'legacy');
            lk = lk + length(source(subdiv).hemisph(1).pinfo);
        end
        pk = cat(2,pk,lk(:)');                 
    end
    
    % Add nearest neighbors to candidate support in subdiv
    %if ~isequal(subdiv,subdiv_sol)
    pk = unique([ico(subdiv).patch(pk).neighbors]);
    %end
end
    
%% Remove Superficial Sources from Final Support
% Find superficial patch numbers (at least 10% of dipoles are < 5mm to inner_skull.surf)
B = [ico(subdiv).patch.mindistout];
ii_s = find(B);

% Form reduced-order lead field (Gtheta) of candidate support (KK)
KK = [];
Gn = [];
Gtheta = [];
currstr = [];

for ii = 1:length(pk)
    % Compute reduced-order lead field of patch pk(ii)
    G_temp = ico(subdiv).patch(pk(ii)).U(:,1:p)*ico(subdiv).patch(pk(ii)).S(1:p,1:p);
    cstr_temp  = repmat(cort_cd*ico(subdiv).patch(pk(ii)).A,p,1);
    
    % Remove superficial patches from support and form reduced-order lead field
    if ~ismember(pk(ii),ii_s)
        KK = cat(2,KK,pk(ii));
        Gn = cat(2,Gn,ico(subdiv).patch(pk(ii)).U(:,1:p));
        Gtheta = cat(2,Gtheta,G_temp);
        currstr = cat(1, currstr, cstr_temp);
    else
        disp(['Excluding superficial patch ' num2str(pk(ii)) ' in ico-' num2str(subdiv) 'p (at least ' num2str(pt_s) '% of dipoles are < 5mm to inner_skull.surf)'])
    end 
end

end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%