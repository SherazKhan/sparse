function [X_greedy, n_greedy, leg_greedy, X_mne, n_mne, leg_mne, X_mne_all, n_mne_all, leg_mne_all] =  ... 
plot_timecourses_megeeg(Xls_plot_grcond, Xls_plot_mnecond, Xls_plot_mnecond_all, time_est, T, sel_rois, leg_mne_all, col, thresh, fullcort, sel_times)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Plot Time Courses and Compute Norms of Resultant Estimates
% Written by Pavitra Krishnaswamy (pavitrak@nmr.mgh.harvard.edu)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% SCSP Scaled Time Series Corresponding to Chosen Regions - Hypersampled Dipole Space
figure, set(gcf,'color','white'); hold on; 
X_greedy = zeros(length(sel_rois),T); n_greedy = zeros(1, length(sel_rois));%initialization 
greedy_col = col;                                                           %initialization
for j = 1:length(sel_rois)
    X_greedy(j,:) = Xls_plot_grcond(j).dip_res(:,1:T);                      %scalarX if Xsim_plot.dip_resovern
    n_greedy(j) = norm(X_greedy(j,sel_times))*10^9;                         %norm of time course [nAm]
    plot(time_est(1:T), smooth(X_greedy(j,:)*10^9),'LineWidth',1.5,'color',greedy_col(j,:)); %plot of time courses
end
xlabel('Time (Milliseconds)','FontSize',14); ylabel('Source Current (nAm)','FontSize',14); 
title('SCSP Dipole Current Estimates','FontSize',14); leg_greedy = sel_rois; legend(leg_greedy,'FontSize',14); set(gca,'FontSize',14);

%% MNE Scaled Time Series - Hypersampled Dipole Space
figure, set(gcf,'color','white'); hold on; 
X_mne = zeros(length(sel_rois),T); n_mne = zeros(1, length(sel_rois));      %initialization
mne_col = col;                                                              %initialization
for j = 1:length(sel_rois)
    X_mne(j,:) = Xls_plot_mnecond(j).dip_res(:,1:T);                        %scalarX if Xsim_plot.dip_resovern
    n_mne(j) = norm(X_mne(j,sel_times))*10^9;                               %norm of time course [nAm]
    plot(time_est(1:T), smooth(X_mne(j,:)*10^9),'LineWidth',1.5,'color',mne_col(j,:)); %plot of time courses
end
xlabel('Time (Milliseconds)','FontSize',14); ylabel('Source Current (nAm)','FontSize',14);
title('MNE Dipole Current Estimates','FontSize',14); leg_mne = sel_rois'; legend(leg_mne,'FontSize',14); set(gca,'FontSize',14); 

%% MNE Scaled Time Series - Hypersampled Dipole Space
X_mne_all = zeros(length(Xls_plot_mnecond_all),T);                          %initialization
n_mne_all = zeros(1,length(Xls_plot_mnecond_all));                          %initialization
for j = 1:length(Xls_plot_mnecond_all)
    X_mne_all(j,:) = Xls_plot_mnecond_all(j).dip_res(:,1:T);                %scalarX if Xsim_plot.dip_resovern
    n_mne_all(j) = norm(X_mne_all(j,sel_times))*10^9;                       %norm of time course [nAm]
end
plot_cols = find(n_mne_all > prctile(n_mne_all,thresh));                    %threshold MNE estimate to 95%ile
figure, set(gcf,'color','white'); hold on;
mneall_col = hsv(length(plot_cols));                                        
for i = 1:size(mneall_col, 1)
    plot(time_est(1:T), smooth(X_mne_all(plot_cols(i),:)*10^9),'LineWidth',1.5,'color',mneall_col(i,:));
end
xlabel('Time (Milliseconds)','FontSize',14); ylabel('Source Current (nAm)','FontSize',14);
title('MNE Dipole Current Estimates','FontSize',14); legend(leg_mne_all(plot_cols),'FontSize',14); set(gca,'FontSize',14); 
if fullcort
    n_mne_all = n_mne_all(plot_cols); leg_mne_all = leg_mne_all(plot_cols);
end

end