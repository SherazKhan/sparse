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
    plot(time*1000, xcest,'color',cortcol(ii,:)); hold on;
    clear ind est
  end

  hold on; 
  for i = 1:length(sim_on)
      line(time(sim_on{i})*1000,35*[1 1],'color','k'); %lines to indicate when simulation is on
  end

  xlabel('Time (Milliseconds)','FontSize',14); ylabel('Source Currents (nAm)','FontSize',14);
  title('Cortical Estimates (ICO-3)','FontSize',14); set(gca,'FontSize',14);
end
end