function pw_dist = pairwise_dist(sras_coord)

%% Initialize
L = length(sras_coord);                                                     % number of regions
avg_coord = zeros(L,3);                                                     % initialize centroids
pw_dist = zeros(L,L);                                                       % pairwise distances

%% Compute Mean Distances B/W Division
for i = 1:L
	for j = 1:L
    	avg_coord(i,:) = mean(sras_coord{i},2);                             % centroid of region 1 [mm]
        avg_coord(j,:) = mean(sras_coord{j},2);                             % centroid of region 2 [mm]
        pw_dist(i,j) = sqrt(sum((avg_coord(i,:)-avg_coord(j,:)).^2,2))/10;  % divide by 10 to go from [mm] to [cm]
	end
end

end