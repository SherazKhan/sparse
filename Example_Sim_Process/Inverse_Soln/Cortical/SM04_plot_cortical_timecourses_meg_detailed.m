%% Initializations
cfwd = mne_read_forward_solution(strcat(datapath,'fwd/', prefix, '-meg-all-fixed-fwd.fif'));         %cort fwd
csrc = mne_read_source_spaces(strcat(datapath,'src/', prefix, '-ico-',num2str(cdiv),'p-src.fif'));   % cort source space
cdiv = 3;                                                                   % cortical division for the key plots 
meas_fname = strcat(datapath,'meas/', prefix,'_vpl_somato_meg_sim',append,'.mat'); % simulated measurement file
load(meas_fname,'Xc_plot','scalarX','stc','t'); alg = 'MNE';                % load inputs into subcortical inverse problem

%% Simulated Time Course Currents: Cortical
figure, set(gcf,'color','white'); hold on; 
num_creg = length(Xc_plot);
Xc = zeros(num_creg,T); cortcol = col(1:num_creg,:);
for i = 1:num_creg
    Xc_series = Xc_plot(i).dip_res;                                         %scalarX if Xc_plot.dip_resovern
    Xc(i,:) = smooth(Xc_series(1:T)*10^9);plot(time_est(1:T), Xc(i,:), 'color',cortcol(i,:),'LineWidth',2); hold on; 
    nc(i) = norm(Xc(i,:));
end
X_sim = Xc; n_sim = nc;                                                     %simulated current time courses
xlabel('Time (Milliseconds)','FontSize',14); ylabel('Source Current (nAm)','FontSize',14);
title('Simulated Activity','FontSize',14); leg_sim = {'cs1','cs2','ppc','is2'}; legend(leg_sim); set(gca,'FontSize',14);


%% Display Result - One Time Series Per Patch, Need Not be As Many Time Series as L
if ploton
for subdiv_plot = 1:3
  figure(3*subdiv_plot), set(gcf,'color','white');
  patchno = patch_est{subdiv_plot};                                         %patch set in ico-subdiv_plot
  cortcol = hsv(length(patchno));                                           %colors for cortical plot
  for ii = 1:length(patchno)
    est = Xest(subdiv_plot).patch{ii};                                      %time course estimates
    [~,ind] = sort(dot(est,est,2),'descend');                               %norms 
    xcest = smooth(est(ind(1),:)*10^9);                                     %time course sorted by norm
    plot(time*1000, xcest,'color',cortcol(ii,:)); hold on;                  %plot time course
    clear ind est
  end
  xlabel('Time (Milliseconds)','FontSize',14); ylabel('Source Currents (nAm)','FontSize',14);
  title('Cortical Estimates (ICO-3)','FontSize',14); set(gca,'FontSize',14);
end
end

%% Assign Labels to Selected Patches - ICO-3
mkdir labels;
cpatch_num = patch_est(cdiv);                                               %final patch set
for i = 1:length(cpatch_num)
    [clblno(i), ~, ~, clr(i)] = patch_to_labloc(csrc,cfwd,cpatch_num(i));   % label # corresponding to patch
    clabel_fname{i} = strcat(num2str(clblno(i),'%06g'),'-',clr(i),'h.label'); %cortical label name
    clabpath = strcat(datapath, 'label/',prefix,'-ico-',num2str(cdiv),'/',clr(i), 'h/'); %cortical label path
    clabsource = strcat(clabpath,clabel_fname{i});        
    copyfile(clabsource,strcat(pwd,'/labels'));
    disp(strcat('Imported cortical label file .... ',clabel_fname{i}));
end
disp('please load labels or submasks into mne-analyze or freeview and confirm on corresponding roi')
rois = {'cs1','cs1','cs2','cs1','is2','cs1','ppc'};                         %rois of selected patches %%%%%

%% Display Result - One Time Series Per Region - Only for ICO-3
figure(subdiv_plot*5), set(gcf,'color','white');
[~, uni_divs, ~] = unique(cpatch_num);                                      %check if this works - else - 1:p:length(...
Xls_plot_grcond = roi_sum(rois_gr, unique_rois, uni_divs, Xls_plot_greedy,'greedy'); % condensed greedy time course estimates by roi

%% Bar Graphs of Current Distribution Across Regions - Only for ICO-3
figure(subdiv_plot*5), set(gcf,'color','white'); 
regions_list = cell(1, length(ico(cdiv).patch)); n_greedy = zeros(1, length(ico(cdiv).patch)); 
for i = 1:length(regions_list)
    regions_list{i} = strcat('patch',num2str(i));                           %regions list
end
for ii = 1:length(cpatch_num)
    n_greedy(cpatch_num) = norm(Xest(subdiv_plot).patch{ii});               %norm
    leg_greedy{ii} = rois;
end
currdistrn_rois(regions_list, n_sim, leg_sim, n_greedy, leg_greedy);        %bar graph for greedy estimates