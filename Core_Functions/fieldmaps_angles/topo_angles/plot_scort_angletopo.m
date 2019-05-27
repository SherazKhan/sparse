function plot_scort_angletopo(prefix, dpath, sdiv, hmgz_fname, pangle, flag)

%% Store Subdiv Indices in sdiv for Mask Access
for k = 1:length(sdiv)
    all_regions{k} = sdiv(k).name;
    if isempty(strfind(all_regions{k},'lhipsurf')) == 0
        ind(k) = 1;
    end
    if isempty(strfind(all_regions{k},'ic')) == 0
        sdiv(k).name = 'lic'; 
        sdiv(end+1) = sdiv(k); 
        sdiv(end).name = 'ric'; sdiv(end).subdiv = 1;
        pangle(end+1) = pangle(k);
    end
end
hipsurf_ind1 = find(ind,1);
list_of_regs = unique(all_regions);
for i = 1:length(list_of_regs)
    ind_sfwd = strncmp(list_of_regs{i}, all_regions, 4);
    ssubdiv = find(ind_sfwd);
    for j = 1:length(ssubdiv)
        sdiv(ssubdiv(j)).subdiv = ssubdiv(j) - ssubdiv(1) + 1;
    end
    clear ind_sfwd
end

%% Assign Angles to Subcortical Volume and Surface Masks
T = 2; newn = zeros(256,256,256,T);                                         %T=2 even though we don't need time for coloring
for i = 1:length(sdiv)
    if isempty(strfind(sdiv(i).name,'hipsurf'))
        % Subcortical Volumes
        maskname = strcat(sdiv(i).name, '_submask',num2str(sdiv(i).subdiv),'.mgz');
    else
        % Hippocampal Surfaces
        maskname = hmgz_fname{i-hipsurf_ind1+1};    
    end
    
    % Assign Angle to Mask
	A = MRIread(strcat(dpath,'mri/svolume_decomps/', maskname));            %subcortical mask
    M = A.vol == 1;                                                         %identify region where mask is high
    f = find(M == 1);                                                       %find indices
    [x, y, z] = ind2sub(size(M),f);                                         %find x y z coordinates in voxel space
    for mm = 1:length(x)
    	for t = 1:T
        	a(t) = pangle(i);                                               %coloring should be different for each subdiv
            newn(x(mm),y(mm),z(mm),t) = a(t);                               %new mask
        end
    end
    clear M f x y z
end

%% Write Full Mask with Angles
B = A; B.vol = newn;                                                        %mask to write
savename = strcat(prefix,'_',flag,'.mgz');                                  %mask name to save
MRIwrite(B,savename)                                                        %write mask in mgz format

end