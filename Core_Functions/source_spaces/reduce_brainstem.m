function [srclocs, ndipoles, volume] = reduce_brainstem(datapath, prefix, orig_srclocs, orig_ndipoles, orig_volume)
% Inputs are Original Values (from Freesurfer) for src locations, # dipoles and volume
% Pavitra Krishnaswamy [pavitrak@nmr.mgh.harvard.edu], Sep 1, 2016
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Load IC
lic_srclocs = importdata(strcat(datapath, '/src/',prefix,'_lic.txt')); 
ric_srclocs = importdata(strcat(datapath, '/src/',prefix,'_ric.txt')); 
ic_srclocs = [lic_srclocs; ric_srclocs];

%% Find Lowest (Most Inferior) Point of IC
lowest_pt  = min(ic_srclocs(:,3));
excl_pts = orig_srclocs(:,3) < lowest_pt;                                   %assess which points in original srclocations are inf to IC
srclocs = orig_srclocs(~excl_pts,:);                                        %consider only midbrain srclocs for future analysis
ndipoles = size(srclocs,1)/3;                                               %update number of dipoles
volume = orig_volume*ndipoles/orig_ndipoles;                                %scale volume by number of dipoles, as 1 dipole = 1 mm3 voxel

end