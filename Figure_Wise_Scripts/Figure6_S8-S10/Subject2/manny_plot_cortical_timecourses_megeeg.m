%% Display Result - One Time Series Per Patch, Need Not be As Many Time Series as L
if ploton
for subdiv_plot = 1:3
  figure(3*subdiv_plot), set(gcf,'color','white');
  patchno = patch_est{subdiv_plot};
  cortcol = hsv(length(patchno));
  for ii = 1:length(patchno)
    est = Xest(subdiv_plot).patch{ii};
    [~,ind] = sort(dot(est,est,2),'descend');
    xcest = smooth(est(ind(1),:)*10^9);
    norm_cort(ii) = norm(xcest); 
    plot(time*1000, xcest,'color',cortcol(ii,:)); hold on;
    clear ind est
  end
  xlabel('Time (Milliseconds)','FontSize',14); ylabel('Source Currents (nAm)','FontSize',14);
  title('Cortical Estimates (ICO-3)','FontSize',14); set(gca,'FontSize',14);
end
end

leg_cort = {'la1','ra1','ra1'};
subc_svds_fname = strcat(datapath, 'fwd/', prefix,'_subc_',meastype,'_plr-prewhite','_SVDs_HFABR.mat'); %subcortical SVD
load(subc_svds_fname,'sregname');                                           %subcortical fwds by subdivisions with regnames
regions_list = [leg_cort, strcat(sregname,'-')];                            %select cortical (based on shortlist from MNE) and all subcort regis
thresh2 = 100;
fullcort = 0;
currdistrn_rois(regions_list, norm_cort, leg_cort, 'scsp', thresh2, fullcort, strcat(sregname,'-')); %bar graph for greedy estimates
