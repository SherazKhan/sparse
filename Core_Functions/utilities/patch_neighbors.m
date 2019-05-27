function ico = patch_neighbors(ico,ii,puse,source,subdiv)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% DEFINE CORTICAL PATCH NEIGHBORHOODS
% Written by Gabriel Obregon-Henao (obregon@nmr.mgh.harvard.edu)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Read triangulations in subdiv
use_tris = source(subdiv).hemisph(1).use_tris;
use_tris2 = source(subdiv).hemisph(2).use_tris;
cols = size(use_tris,2);

%% Find triangles containing iith vertex in left and right hemispheres (tri and tri2, respectively)
tri = reshape(use_tris',1,size(use_tris,1)*cols);
tri2 = reshape(use_tris2',1,size(use_tris2,1)*cols);
tri = find(ismember(tri,source(subdiv).hemisph(1).vertno(ii)));
tri = ceil(tri/cols);
tri2 = find(ismember(tri2,source(subdiv).hemisph(2).vertno(ii)));
tri2 = ceil(tri2/cols);

%% Find unique vertices conforming these triangles (neighborhood vertices)
vert = unique(use_tris(tri,:));
vert2 = unique(use_tris2(tri2,:));

%% Find corresponding patch numbers
[~,~,ico(subdiv).patch(ii).neighbors] = intersect(vert,source(subdiv).hemisph(1).vertno,'legacy');
[~,~,lk2] = intersect(vert2,source(subdiv).hemisph(2).vertno,'legacy');
ico(subdiv).patch(ii+puse).neighbors = lk2 + puse;
end