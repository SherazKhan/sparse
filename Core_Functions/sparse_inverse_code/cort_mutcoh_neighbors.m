function u = cort_mutcoh_neighbors(svds_fname,p,pt_s,subdiv_stop)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% COMPUTE AVERAGE WORST-CASE MUTUAL COHERENCE B/W ANY TWO RANDOMLY-SELECTED PATCHES
% Written by Gabriel Obregon-Henao (obregon@nmr.mgh.harvard.edu)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

load(svds_fname,'ico');

for subdiv = 1:min(subdiv_stop,4)
    disp(['Computing mutual coherences in ico-' num2str(subdiv) 'p source space'])
    
    %% Merge 1:p(subdiv) modes of all disjoint patches
    Gn = [];   
    for ii = 1:length(ico(subdiv).patch)
        Gn = cat(2,Gn,ico(subdiv).patch(ii).U(:,1:p(subdiv)));
    end
    
    %% Compute mutual coherence matrix (Gc)
    Gc = abs(Gn'*Gn);
    
    %% Find superficial patch numbers (at least <pt_s> percent of dipoles are < 5mm to inner_skull.surf)
    B = [ico(subdiv).patch.mindistout];
    ii_s = find(B);
    
    %% Compute number of patches in both hemispheres
    KK = length(ico(subdiv).patch);
    
    %% Check that iith patch is not superficial and compute maximum mutual coherence with each non-superficial neighbor
    for ii = 1:KK
        max_corr = [];
        if ~ismember(ii,ii_s)
            % Load nearest neighbors
            neigh = setdiff(ico(subdiv).patch(ii).neighbors,ii);
            % Compute maximum mutual coherence between iith patch and jjth neighbor
            for jj = 1:length(neigh)
                if ~ismember(neigh(jj),ii_s)
                    max_corr = cat(2,max_corr,max(max(Gc(ii*p(subdiv)-(p(subdiv)-1):ii*p(subdiv),neigh(jj)*p(subdiv)-(p(subdiv)-1):neigh(jj)*p(subdiv)))));
                end
            end
            % Compute average worst-case mutual coherence between iith patch and all of its nearest neighbors
            mean_corr(ii) = mean(max_corr);
        else
            mean_corr(ii) = 0;
            disp(['Excluding superficial patch ' num2str(ii) ' (at least ' num2str(pt_s) '% of dipoles are < 5mm to inner_skull.surf)'])
        end
    end
    
    %% Compute average worst-case mutual coherence across all neighborhoods
    u(subdiv) = round(mean(mean_corr(mean_corr>0))*10^2)/10^2;   
end
end