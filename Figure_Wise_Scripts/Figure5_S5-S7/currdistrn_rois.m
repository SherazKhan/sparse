function currdistrn_rois(regions_list, n_sim, leg_sim, plot1_leg, n_est, leg_est, plot2_leg, ...
         cort_scort_flag, simon)

num_divs = 3;                                                               % Set up Number of Chars to Compare

%% Quantify RMS Values By Division: Simulation
for m = 1:length(regions_list)
    reg_sim = find(strncmp(leg_sim, regions_list{m}, num_divs));            %all instances of mth region in current source space
    if ~isempty(reg_sim)
        norms_sim = n_sim(reg_sim);                                         %norms of selected regions
        totnorm_sim(m) = sqrt(sum(norms_sim.^2));                           %total norm across selected regions
    else
        totnorm_sim(m) = 0;                                                 %no representation of this region in sim
    end
end

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
ind_sel = 1:length(regions_list);
disp(regions_list(ind_sel));
[~,tmp,~] = unique(regions_list(ind_sel)); ind_sel = ind_sel(tmp); 

%% Rearrange Labels by Cortical and Subcortical
csflag_sel = cort_scort_flag(ind_sel);
sregions_list = regions_list(ind_sel);      rearrange_regions_list = sregions_list([find(csflag_sel) find(~csflag_sel)]);
stotnorm = totnorm_sim(ind_sel);            rearrange_totnorm_sim = stotnorm([find(csflag_sel) find(~csflag_sel)]);
etotnorm = totnorm_est(ind_sel);            rearrange_totnorm_est = etotnorm([find(csflag_sel) find(~csflag_sel)]);

%% Plot Simulation Only
if simon
figure,  set(gcf,'color','white');
bar(rearrange_totnorm_sim);
sim = gca;                       sim.XTick = 1:1:length(ind_sel);
sim.XTickLabel = rearrange_regions_list;
sim.XTickLabelRotation = 45;     sim.YTickLabelRotation = 45;
set(sim,'FontSize',14);
xlabel('Region','FontSize',14); ylabel('RMS Dipole Amplitudes (nAm)', 'FontSize',14)
title([plot1_leg, ' Across Regions'],'FontSize',14); set(gca,'FontSize',14);
legend(plot1_leg)
end

%% Plot Figure
norms_toplot = [rearrange_totnorm_sim' rearrange_totnorm_est'];
figure, set(gcf,'color','white');
bar(norms_toplot);
hb = gca;                       hb.XTick = 1:1:length(ind_sel);
hb.XTickLabel = rearrange_regions_list;
hb.XTickLabelRotation = 45;     hb.YTickLabelRotation = 45;
set(hb,'FontSize',14);
xlabel('Region','FontSize',14); ylabel('RMS Dipole Amplitudes (nAm)', 'FontSize',14)
title(['Compare ', plot1_leg, ' vs. ', plot2_leg, ' Across Regions'],'FontSize',14); set(gca,'FontSize',14);
legend(plot1_leg, plot2_leg)

end