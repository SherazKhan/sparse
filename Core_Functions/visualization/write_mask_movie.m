function write_mask_movie(trimP, mask, tstep, Tst, Ten, Xls_plot, mask_offset, mask_mean, ...
         aseg, fullcort, prefix, meastype, scalarX)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Write Dipole Current Estimate Time Series Across Regions into 1 MRI Mask 
% Can Load the Resulting Mask into Freeview to see Movie of Estimates
% Written by Pavitra Krishnaswamy (pavitrak@nmr.mgh.harvard.edu)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Initialize
t = Tst:tstep:Ten; l = length(t);                                           % time points in estimates
mask_movie = zeros([size(aseg.vol),l]);                                     % initialize mask
offs = min(mask_offset);                                                    % offset of Xls_dip across regions
m = mean(mask_mean);                                                        % mean of Xls_dip across regions

%% Generate Mask Movie
for i = 1:length(trimP)
    Mreg = mask{i};                                                         % focus on mask for region trimP(i)
    f = find(Mreg==1);                                                      % find indices
    [x, y, z] = ind2sub(size(Mreg),f);                                      % find x y z coordinates in voxel space
    for mm = 1:length(x)                                                    % index of dipole (point) in static mask
        for tt = 1:l                                                        % time indices for dipole of interest
            mask_amp = Xls_plot(i).dip_res(t(tt));                          % amplitude of mask - - scalarX if dip_resovern
            %mask_amp = Xls_plot(i).dip_res(t(tt))/m;                       % scaled amplitude of mask - scalarX if dip_resovern
            %mask_amp = Xls_plot(i).dip_res(t(tt))-offs)/m;                 % scaled and shifted amplitude of mask - scalarX if dip_resovern
            mask_movie(x(mm),y(mm),z(mm),tt) = mask_amp;                    % new mask
        end
    end
    clear Mreg f x y z mm;                                                  % put into appropriate format
end

%% Save Mask Movie
B = aseg; B.vol = mask_movie;
if fullcort == 0
    savename = strcat(prefix, '_scsp_',meastype, '_estimates_mask.mgz');    % scsp = subcortical 2 stage subspace pursuit
else
    savename = strcat(prefix, '_sp_',meastype,'_estimates_mask.mgz');       % sp = subspace pursuit
end
MRIwrite(B,savename);                                                       % write MRI mask for scsp estimates
disp('Wrote Mask Movie to Display Greedy Estimates'); 

end