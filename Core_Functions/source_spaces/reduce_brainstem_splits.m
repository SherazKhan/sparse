function [srclocs, ndipoles, volume, pick_inds] = reduce_brainstem_splits(...
         orig_srclocs, orig_ndipoles, orig_volume, num_bs, savename, ploton)
% Inputs are Original Values (from Freesurfer) for src locations, # dipoles and volume
% Pavitra Krishnaswamy [pavitrak@nmr.mgh.harvard.edu], May 25, 2016
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Split into Midbrain, Pons and Medulla
[split_ind, ~] = split_subdiv(orig_srclocs, num_bs, savename, ploton);      %split brainstem into 4 pieces

%% Find Which Split is Midbrain
for j = 1:length(split_ind)
    Z{j} = orig_srclocs(split_ind{j},3);                                    %coordinates in inf-sup direction 
    m(j) = max(Z{j});                                                       %figure out most sup. point of this split
end
[~, pick_split] = max(m);                                                   %pick split which has the most sup. point - i.e. midbrain
pick_inds = split_ind{pick_split};                                          %srcspace indices to pick
srclocs = orig_srclocs(pick_inds,:);                                        %consider only midbrain srclocs for future analysis
ndipoles = size(srclocs,1)/3;                                               %update number of dipoles
volume = orig_volume*ndipoles/orig_ndipoles;                                %scale volume by number of dipoles, as 1 dipole = 1 mm3 voxel

end