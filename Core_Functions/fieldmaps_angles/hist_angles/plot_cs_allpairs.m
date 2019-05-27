function plot_cs_allpairs(angcs_ov, sdiv, cpatch, saveon, chtype, choice_fn)

%% Collate Principal Angle Distributions
all_ang_same = []; all_ang_opp = [];
for i = 1:size(angcs_ov,1)                                                  %subcortical index
    for j = 1:size(angcs_ov,2)                                              %cortical index
        disp(i); disp(j);
        
        % Ipsilateral or Contralateral
        subc_hemi = sdiv(i).hemi;                                           %subcortical hemisphere
        cort_hemi = cpatch(j).hemi;                                         %cortical hemisphere
        if strcmp(subc_hemi, cort_hemi)
            same_hemisph(i,j) = 1;                                          %cort, scort in same hemisphere
            all_ang_same = [all_ang_same angcs_ov{i,j}(:)'];                %all the angles for ipsilateral
        elseif ~strcmp(subc_hemi, cort_hemi)
            same_hemisph(i,j) = 0;                                          %cort, scort in opposite hemispheres
            all_ang_opp = [all_ang_opp angcs_ov{i,j}(:)'];                  %all the angles for contralateral
        end
        
        % Angle Distributions
        meanov_ang(i,j) = mean(angcs_ov{i,j});                              %mean principal angle for all representative modes
        if ~isempty(angcs_ov{i,j})
        	minov_ang(i,j) = min(angcs_ov{i,j});                            %min principal angle for all representative modes
            maxov_ang(i,j) = max(angcs_ov{i,j});                            %max principal angle for all representative modes
        else
            minov_ang(i,j) = nan;                                           %empty
            maxov_ang(i,j) = nan;                                           %empty
        end
    end
end

%% Plot Histograms
if ~isempty(minov_ang)
for choose_fn = choice_fn:1:choice_fn                           
    % Collate Min, Mean, Max or All Angles
    switch choose_fn 
        case 0
            asame = minov_ang(same_hemisph==1);     asame = asame(~isnan(asame));  
            aopp = minov_ang(same_hemisph==0);      aopp = aopp(~isnan(aopp));
            savename = strcat('cs_',chtype,'_min_principal_angles');
        case 1
            asame = meanov_ang(same_hemisph==1);    asame = asame(~isnan(asame));
            aopp = meanov_ang(same_hemisph==0);     aopp = aopp(~isnan(aopp));
            savename = strcat('cs_',chtype,'_mean_principal_angles');
        case 2
            asame = maxov_ang(same_hemisph==1);     asame = asame(~isnan(asame));
            aopp = maxov_ang(same_hemisph==0);      aopp = aopp(~isnan(aopp));
            savename = strcat('cs_',chtype,'_max_principal_angles');
        case 3
            asame = all_ang_same;       
            aopp = all_ang_opp;
            savename = strcat('cs_',chtype,'_all_principal_angles');
        case 4
            a_all = [all_ang_same all_ang_opp];
            savename = strcat('cs_',chtype,'_all_principal_angles_nohemis');
    end
    
    % Plot Histograms of Collated Angles
    figure, set(gcf,'color','white');
    if choose_fn < 4
      subplot(2,1,1), [~,Xs] = hist(asame);
      histnorm(asame,Xs); %hist(asame); 
      xlabel('Principal Angle (Degrees)','FontSize',14); ylabel('Normalized Histogram','FontSize',14); 
      legend({'Ipsilateral'},'FontSize',14); xlim([0 90]); set(gca,'FontSize',14);    
      subplot(2,1,2), [~,Xo] = hist(aopp);
      histnorm(aopp,Xo); %hist(aopp); 
      xlabel('Principal Angle (Degrees)','FontSize',14); ylabel('Normalized Histogram','FontSize',14); 
      legend({'Contralateral'},'FontSize',14); xlim([0 90]); set(gca,'FontSize',14); 
      clear choose_ang asame aopp Xs Xo
    else
      [~,Xs] = hist(a_all);
      histnorm(a_all,Xs); %hist(a_all);
      xlabel('Principal Angle (Degrees)','FontSize',14); ylabel('Normalized Histogram','FontSize',14); 
      xlim([0 90]); set(gca,'FontSize',14);        
    end
    
    % Save Histograms of Collated Angles
    if saveon
        print(gcf,'-depsc',savename); 
    end
end
end

end