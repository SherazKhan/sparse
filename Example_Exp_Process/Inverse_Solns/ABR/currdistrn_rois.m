function currdistrn_rois(regions_list, n_est, leg_est, plot2_leg, thresh, fullcort, regname)

num_divs = 3;                                                               % Set up Number of Chars to Compare

%% Quantify RMS Values By Division: Estimation
for m = 1:length(regions_list)
    reg_est = find(strncmp(leg_est, regions_list{m},num_divs));             %all instances of mth region in current source space
    if ~isempty(reg_est)
        norms_est = n_est(reg_est);                                         %norms of selected regions
        totnorm_est(m) = sqrt(sum(norms_est.^2));                           %total norm across selected regions
    else
        totnorm_est(m) = 0;                                                 %no representation of this region in sim
    end
end

%% Select Labels to Plot
if fullcort
    ind_sel1 = find(totnorm_est > prctile(totnorm_est,thresh));
    ind_sel2 = find(ismember(regions_list, regname));
    ind_sel = unique([ind_sel1, ind_sel2]);
else
    ind_sel = 1:length(regions_list);
end
disp(regions_list(ind_sel));
[~,tmp,~] = unique(regions_list(ind_sel)); ind_sel = ind_sel(tmp);

%% Plot Figure
norms_toplot = totnorm_est';
figure, set(gcf,'color','white');
bar(norms_toplot(ind_sel,:));
hb = gca;                       hb.XTick = 1:1:length(ind_sel);
hb.XTickLabel = regions_list(ind_sel);
hb.XTickLabelRotation = 45;     hb.YTickLabelRotation = 45;
set(hb,'FontSize',14);
xlabel('Region','FontSize',14); ylabel('RMS Dipole Amplitudes (nAm)', 'FontSize',14)
title([plot2_leg, ' Across Regions'],'FontSize',14); set(gca,'FontSize',14);
legend(plot2_leg)

end