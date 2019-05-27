function [xs_scaled, Gs_unscaled, xfit_scaled, Gfit_unscaled, lkmax, ys] = ...
         fit_scort_to_sparsecort(sdiv, subset_inds, SNR_fit, s_currstr, c_currstr)

%% Simulate Activity in Most Significant Eigenmode
smodes = subset_inds.smodes;                                                % Mode of Interest - should be 1!
Us = sdiv.U(:,smodes);                                                      % Fwd Mode of Interest (Whitened)
Ss = sdiv.S(smodes,smodes);                                                 % Fwd Spectral Norm of Interest (Whitened)
Gs = s_currstr*Us*Ss;                                                       % Prewhitened Subcortical Fwd: Reduced Dim + Scaled by Currstr [Am]
if length(smodes) == 1
    xs = 1;                                                                 % Source Currents: Mode Space - 1 unit of currstr of this division (No Units)
else
    xs = 1./sqrt(diag(Ss)*sum(diag(Ss)));                                   % Source Currents: Mode Space - xs propto SingVals (No Units)
end
ys = Gs*xs;                                                                 % Data for Field Map : Prewhitened (No Units)
N = size(ys,1);                                                             % # Sensors
xs_scaled = xs*s_currstr;                                                   % scaled [Am]
Gs_unscaled = Us*Ss;                                                        % for output if use xx_scaled

%% Resultant Amplitude of Simulated Activity in Subcortical Region
Vs = sdiv.V(:,smodes);                                                      % Right eigenvector - num_dipoles X 1
Xs_plot.mode = xs*s_currstr;                                                % Mode Amplitudes [nAm]
Xs_plot.dip_allorient =  Vs*Xs_plot.mode;                                   % Dipole Amplitudes [Am]
Xs_plot.dip_or(:,:,1) = Xs_plot.dip_allorient(1:3:end,:);                   % Orientation 1
Xs_plot.dip_or(:,:,2) = Xs_plot.dip_allorient(2:3:end,:);                   % Orientation 2
Xs_plot.dip_or(:,:,3) = Xs_plot.dip_allorient(3:3:end,:);                   % Orientation 3
[Xs_plot.dip_res, Xs_plot.dip_resor] = vec_sum(Xs_plot.dip_or);             % Net simulated resultant activity for region [Am]
Xs_source = Xs_plot.dip_res;                                                % Resultant Source Current [Am]

%% Collate Forward Solutions of All Cortical Patches in Sparse Circuit
Uc = subset_inds.cfwd_emode;                                                % Fwd Modes of Interest (Whitened)
Sc = diag(subset_inds.cfwd_eval);                                           % Fwd Spectral Norm of Interest (Whitened)
Gc = c_currstr*Uc*Sc;                                                       % Reduced Cortical Forward Solution (not scaled by currstr)
%the above assumes none of the patches in the circuit are medial wall patches, and hence dont need to be excluded
Gfit = Gc; lkmax = size(Gfit,2);                                            % Number of cortical modes to consider for fit
Gfit_unscaled = Uc*Sc;                                                      % for output if use xfit_scaled

%% Fit Activity to All Cortical Patches
xfit = source_estimates('MNE', eye(N), Gfit, N, ys, 1, size(ys,2), SNR_fit, [], ones(size(Gfit,2))); %Fitted Mode Amplitudes (no Units)
xfit_scaled = xfit*c_currstr;                                               % Fitted Mode Amplitudes [nAm]
%Xc_plotmri = srcspace_activity_mri(prefix, L, xfit_scaled, cortp, cfwd, ico_num, source, ico, mri_path, patch_exc, lkmax, fittype, Xs_source);
%scaled_patchres = [Xc_plotmri(:).dip_res_sc];                               % Mode to Dip, account orient, compute patch-wise resultants, scale by Xs_source
% figure, set(gcf,'color','white'); plot(scaled_patchres); 

end