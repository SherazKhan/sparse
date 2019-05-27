function SM04_simplotfinal(overlap)
% PLOT SIMULATED ACTIVITY and MEASUREMENT

%% Load Data
datapath = '/autofs/eris/purdongp/users/pavitrak/sourceloc_data/SEP/SM04/';
measname = strcat(datapath, 'meas/SM04_vpl_somato_meg_sim',overlap,'.mat'); 
load(measname); 
time = time*1000;                                                           %time in milliseconds
cortcol = [0 0 0; 1 0 0; 0 1 0; 0.5 0.5 1];                                 %colors for cortical regions
stc_rewrite = 0;                                                            %rewrite stc file
num_creg = length(Xc_plot);                                                 %number of cortical regions

%% Source Currents (Resultant)
figure(1), set(gcf,'color','white');
% Subcortical
Xvpl_series = Xvpl_plot.dip_res;                                            %scalarX if Xvpl_plot.dip_resovern
[~,as] = sort(dot(Xvpl_series,Xvpl_series,2),'descend');
Xs = smooth(Xvpl_series(as(1),:)*10^9); plot(time, Xs,'m'); hold on; 
for i = 1:length(sim_on)
    line(time(sim_on{i}),65*[1 1],'color','k');
end
% Cortical
Xc = zeros(num_creg,t); 
for i = 1:num_creg
    Xc_series = Xc_plot(i).dip_res;                                         %scalarX if Xc_plot.dip_resovern
    [~,ac] = sort(dot(Xc_series,Xc_series,2),'descend');
    Xc(i,:) = smooth(Xc_series(ac(1),:)*10^9);
    plot(time, Xc(i,:), 'color',cortcol(i,:)); hold on; 
end
% Plot Settings
set(gca,'FontSize',14);
xlabel('Time (Milliseconds)','FontSize',14); ylabel('Source Current (nAm)','FontSize',14);
title('Simulated Activity','FontSize',14); legend('VPL Thalamus','cS1','cS2','PPC','iS2');
print -depsc -r300 -painters sim_source_curr

%% Noise Free Measurement
figure(2), set(gcf,'color','white'); 
subplot(2,1,1), plot(time, Ynoisefree(magnetom,:)*10^15); 
set(gca,'FontSize',14);
title('Noiseless Measurement','FontSize',14);       
xlabel('Time (Milliseconds)','FontSize',14); ylabel('Magnetometer (fT)','FontSize',14);
subplot(2,1,2), plot(time, Ynoisefree(gradiom,:)*10^13); 
set(gca,'FontSize',14);
xlabel('Time (Milliseconds)','FontSize',14); ylabel('Gradiometer (fT/cm)','FontSize',14);
print -depsc -r300 -painters noiseless_meas

%% Noisy Measurement
figure(3), set(gcf,'color','white');    
subplot(2,1,1), plot(time, Y(magnetom,:)*10^15); 
title('Noisy Measurement','FontSize',14);
set(gca,'FontSize',14);
xlabel('Time (Milliseconds)','FontSize',14); ylabel('Magnetometer (fT)','Fontsize',14);
subplot(2,1,2), plot(time, Y(gradiom,:)*10^13);
set(gca,'FontSize',14);
xlabel('Time (Milliseconds)','FontSize',14); ylabel('Gradiometer (fT/cm)','FontSize',14);
print -depsc -r300 -painters noisy_meas

%% View Simulated STCs (mne_analyze) to Ensure Correct Regions and Consistent Scales
if stc_rewrite
    X_lh = cat(1,X(1:3).data); M = size(X_lh,1);
    stc(1).data = Xhyper_plot.res(1:M,:);                                   %scalarX if Xhyper_plot.resovern
    stc(2).data = Xhyper_plot.res(M+1:end,:);                               %scalarX if Xhyper_plot.resovern
    mne_write_stc_file1(strcat(measname,'-lh.stc'),stc(1))
    mne_write_stc_file1(strcat(measname,'-rh.stc'),stc(2))
end

%% Covariances
figure(4), set(gcf,'color','white');
imagesc(C0(magnetom,magnetom));
title('Magnetometer Covariance','FontSize',14); 
xlabel('Magnetometer Index','FontSize',14); ylabel('Magnetometer Index','FontSize',14); 
colorbar; set(gca,'FontSize',14);

figure(5), set(gcf,'color','white');
imagesc(C0(gradiom,gradiom));
title('Gradiometer Covariance','FontSize',14); 
xlabel('Gradiometer Index','FontSize',14); ylabel('Gradiometer Index','FontSize',14); 
colorbar; set(gca,'FontSize',14);

%% Compute RMS Values for Bar Graph of Simulation Norms Across Regions
for i = 1:num_creg
    sim_Xest = cat(1,Xc_plot(i).dip_res(:,1:t));
    cort_sim_norm(i) = sqrt(dot(sim_Xest,sim_Xest,2));
end
clear sim_Xest
sim_Xest = cat(1,Xvpl_plot.dip_res(:,1:t));
vpl_sim_norm = sqrt(dot(sim_Xest,sim_Xest,2));
clear sim_Xest

%% Plot Bar Graphs
regions_list = {'cortical','lp','lc','lt','llgn','lmgn','ic',...
                'rp','rc','rt','rlgn','rmgn','lhipsurf','rhipsurf'};
sim_vals = zeros(1,length(regions_list)); %t = 40;
sim_vals(1) = mean(cort_sim_norm);
sim_vals(4) = mean(vpl_sim_norm);
figure(6), set(gcf,'color','white');
norms_toplot = [sim_vals'];
bar(norms_toplot*10^9); set(gca,'XTick',1:14,'XTickLabel',regions_list,'FontSize',14);
xlabel('Region','FontSize',14); ylabel('RMS Dipole Amplitudes (nAm)', 'FontSize',14)
title('Simulated Norms Across Regions','FontSize',14); set(gca,'FontSize',14);

%% Bar Graph Computing Norms of Simulation
regions_list2 = {'cS1', 'cS2', 'PPC', 'iS2','lp','lc','lt','llgn','lmgn','ic',...
                'rp','rc','rt','rlgn','rmgn','lhipsurf','rhipsurf'};
sim_vals2 = zeros(1,length(regions_list2)); %t = 40;
sim_vals2(1:4) = cort_sim_norm;
sim_vals2(7) = mean(vpl_sim_norm);
figure(7), set(gcf,'color','white');
norms_toplot2 = [sim_vals2'];
bar(norms_toplot2*10^9); set(gca,'XTick',1:17,'XTickLabel',regions_list2,'FontSize',14);
xlabel('Region','FontSize',14); ylabel('RMS Dipole Amplitudes (nAm)', 'FontSize',14)
title('Simulated Norms Across Regions','FontSize',14); set(gca,'FontSize',14);
end