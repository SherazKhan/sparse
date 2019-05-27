function u = scort_mutcoh_neighbors(svds_fname, varargin)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% COMPUTE AVERAGE WORST-CASE MUTUAL COHERENCE B/W TWO RANDOMLY-SELECTED SUBDIVISIONS
% Written by Gabriel Obregon-Henao (obregon@nmr.mgh.harvard.edu) 
% Written by Pavitra Krishnaswamy (pavitrak@nmr.mgh.harvard.edu)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

load(svds_fname,'sdiv_fwd','sregname');
for i = 1:length(sdiv_fwd)
    all_regions{i} = sdiv_fwd(i).reg_name;
end
if nargin > 1
    ind_sfwd = zeros(1, length(sdiv_fwd));
    excl_regnames = zeros(1, length(sregname));
    names_regs = varargin{1};
    for i = 1:length(names_regs)
      ind_sfwd = ind_sfwd + strncmp(names_regs{i}, all_regions, 4); 	   %find indices of sdiv_fwd that match varargin{i}
      excl_regnames = excl_regnames + strncmp(names_regs{i}, sregname, 4); %find indices of sregname that match varargin{i}
    end
    sdiv_fwd = sdiv_fwd(~ind_sfwd);
    all_regions = all_regions(~ind_sfwd);
    sregname = sregname(~excl_regnames);
end

%% Collate All Region Names
num_regions = length(sregname);
mean_corr = zeros(num_regions, 1);
for i = 1:length(sdiv_fwd)
    all_regions{i} = sdiv_fwd(i).reg_name;
end

%% Compute Max. Mutual Coherence with Each Subdivision within Anatomic Regions
for ii = 1:num_regions
    % Divisions within Current Anatomic Region (Defines a Neighborhood)
    disp(['Computing mutual coherences in neighborhood of ...', sregname{ii}]);
    subreg_ind = find(ismember(all_regions, sregname{ii}));
    num_subregions = length(subreg_ind);
    max_corr = zeros(num_subregions);
    
    for jj = 1:num_subregions
        % Read in modes of forward matrix for jjth subdiv
        sdiv_jj = sdiv_fwd(subreg_ind(jj));
        Gjj = sdiv_jj.whiteUs(:,1:sdiv_jj.numwhite_modes);
        
        % Compute maximum (worst case) mutual coherence between jjth subdiv & kkth neighbor subdivs
        for kk = jj+1:num_subregions
            sdiv_kk = sdiv_fwd(subreg_ind(kk));
            Gkk = sdiv_kk.whiteUs(:,1:sdiv_kk.numwhite_modes);
            max_corr(jj,kk) = max(max(Gjj'*Gkk));
            max_corr(kk,jj) = max_corr(jj,kk);
        end
    end
     % Compute average worst-case mutual coherence b/w all neighbor subdivs within iith region
     mean_corr(ii) = mean(mean(max_corr(max_corr>0)));
end

%% Compute average worst-case mutual coherence across all neighborhoods    
u = round(mean(mean_corr)*100)/100;

end