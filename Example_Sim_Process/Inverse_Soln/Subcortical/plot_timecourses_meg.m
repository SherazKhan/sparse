function [X_sim, n_sim, leg_sim, cflag_sim, X_greedy, n_greedy, leg_greedy, cflag_greedy, X_mne, n_mne, leg_mne, cflag_mne, ...
X_mne_all, n_mne_all, leg_mne_all, cflag_mne_all] =  plot_timecourses_meg(Xvpl_plot, Xc_plot, Xls_plot_grcond, Xls_plot_mnecond, ...
Xls_plot_mnecond_all, time_est, T, sel_rois, sel_flag, leg_all, cflag_all, col, thresh, fullcort, sel_times, sim_on)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Plot Time Courses and Compute Norms of Resultant Estimates
% Written by Pavitra Krishnaswamy (pavitrak@nmr.mgh.harvard.edu)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Simulated Time Course Currents
figure, set(gcf,'color','white'); hold on; 
% Cortical
num_creg = length(Xc_plot);
Xc = zeros(num_creg,T); cortcol = col(1:num_creg,:);
for i = 1:num_creg
    Xc_series = Xc_plot(i).dip_res;                                         %scalarX if Xc_plot.dip_resovern
    Xc(i,:) = smooth(Xc_series(1:T));plot(time_est(1:T), Xc(i,:)*10^9, 'color',cortcol(i,:),'LineWidth',2); hold on; 
    nc(i) = norm(Xc(i,sel_times));
end
% Subcortical
Xvpl_series = Xvpl_plot.dip_res;                                            %scalarX if Xvpl_plot.dip_resovern
Xs = smooth(Xvpl_series(:,1:T)); hold on; plot(time_est(1:T), Xs*10^9,'m','LineWidth',4,'color','m');
hold on;
for i = 1:length(sim_on)
    line(time_est(sim_on{i}),65*[1 1],'color','k');			    %lines to indicate when simulation is on
end
ns = norm(Xs(sel_times));                                                   %norm of time course [nAm]
X_sim = [Xc; Xs']; n_sim = [nc ns]*10^9;                                    %simulated current time courses
xlabel('Time (Milliseconds)','FontSize',14); ylabel('Source Current (nAm)','FontSize',14);
title('Simulated Activity','FontSize',14); leg_sim = {'cs1','cs2','ppc','is2','vpl'}; legend(leg_sim); set(gca,'FontSize',14);
cflag_sim = [1 1 1 1 0];

%% SCSP Scaled Time Series Corresponding to Chosen Regions - Hypersampled Dipole Space
figure, set(gcf,'color','white'); hold on; 
X_greedy = zeros(length(sel_rois),T); n_greedy = zeros(1, length(sel_rois));%initialization 
greedy_col = col;                                                           %initialization
for j = 1:length(sel_rois)
    X_greedy(j,:) = Xls_plot_grcond(j).dip_res(:,1:T);                      %scalarX if Xsim_plot.dip_resovern
    n_greedy(j) = norm(X_greedy(j,sel_times))*10^9;                         %norm of time course [nAm]
    plot(time_est(1:T), smooth(X_greedy(j,:)*10^9),'LineWidth',1.5,'color',greedy_col(j,:)); %plot of time courses
end
hold on;
for i = 1:length(sim_on)
    line(time_est(sim_on{i}),65*[1 1],'color','k');			    %lines to indicate when simulation is on
end
xlabel('Time (Milliseconds)','FontSize',14); ylabel('Source Current (nAm)','FontSize',14); 
title('SCSP Dipole Current Estimates','FontSize',14); leg_greedy = sel_rois; legend(leg_greedy,'FontSize',14); set(gca,'FontSize',14);
cflag_greedy = sel_flag;

%% MNE Scaled Time Series - Hypersampled Dipole Space
figure, set(gcf,'color','white'); hold on; 
X_mne = zeros(length(sel_rois),T); n_mne = zeros(1, length(sel_rois));      %initialization
mne_col = col;                                                              %initialization
for j = 1:length(sel_rois)
    X_mne(j,:) = Xls_plot_mnecond(j).dip_res(:,1:T);                        %scalarX if Xsim_plot.dip_resovern
    n_mne(j) = norm(X_mne(j,sel_times))*10^9;                               %norm of time course [nAm]
    plot(time_est(1:T), smooth(X_mne(j,:)*10^9),'LineWidth',1.5,'color',mne_col(j,:)); %plot of time courses
end
hold on;
for i = 1:length(sim_on)
    line(time_est(sim_on{i}),15*[1 1],'color','k');			    %lines to indicate when simulation is on
end
xlabel('Time (Milliseconds)','FontSize',14); ylabel('Source Current (nAm)','FontSize',14);
title('MNE Dipole Current Estimates','FontSize',14); leg_mne = sel_rois'; legend(leg_mne,'FontSize',14); set(gca,'FontSize',14); 
cflag_mne = sel_flag;

%% MNE Scaled Time Series - Hypersampled Dipole Space
X_mne_all = zeros(length(Xls_plot_mnecond_all),T);                          %initialization
n_mne_all = zeros(1,length(Xls_plot_mnecond_all));                          %initializati
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
hold on;
for i = 1:length(sim_on)
    line(time_est(sim_on{i}),5*[1 1],'color','k');			    %lines to indicate when simulation is on
end
xlabel('Time (Milliseconds)','FontSize',14); ylabel('Source Current (nAm)','FontSize',14);
title('MNE Dipole Current Estimates','FontSize',14); legend(leg_all(plot_cols),'FontSize',14); set(gca,'FontSize',14); 
if fullcort
    n_mne_all = n_mne_all(plot_cols); 
    leg_mne_all = leg_all(plot_cols);
    cflag_mne_all = cflag_all(plot_cols); 
else
    leg_mne_all = leg_all;
    cflag_mne_all = cflag_all;
end

end