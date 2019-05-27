function [Dind, D] = split_subdiv(A, num_divs, savefigname, ploton)
% Written by Pavitra Krishnaswamy and Gabriel Obregon

%% Divide source space into k clusters using K-means Clustering Algorithm
if num_divs > 1
    num_divs = round(num_divs); 
else
    num_divs = 1;
end
src_verts = A(1:3:end,:);
IDX = kmeans(src_verts,num_divs,'Replicates',5,'MaxIter',250,'display','final');
count = 1;
for i = 1:length(IDX)
    IDX_OR(count:count+2) = repmat(IDX(i),1,3);
    count = count + 3;
end

%% Indices of Each Section - Easiest and most Accurate
for i = 1:num_divs
    Dind{i} = find(IDX_OR == i); 
    D{i} = A(Dind{i},:);
end

if ploton
%% Plot Sources to Visualize Clusters
col = hsv(num_divs);
figure, set(gcf,'color','white'); hold all;
for i = 1:num_divs
    scatter3(A(Dind{i},1),A(Dind{i},2),A(Dind{i},3), 10, col(i,:))
end
savefig(strcat(savefigname, '_srcspace_new'));

%% Compute and Plot # Sources per Cluster
for i = 1:num_divs
    num_dipl(i) = sum(ismember(IDX,i));
end
figure, set(gcf,'color','white'); hist(num_dipl);
savefig(strcat(savefigname, '_volumedist_new'));
end

end