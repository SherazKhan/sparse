function plot_cc_allpairs(angcc_ov, cpatch, saveon, chtype, choice_fn)

%% Collate Principal Angle Distributions
all_ang_same = []; all_ang_opp = [];
minov_ang = []; meanov_ang = []; maxov_ang = [];
for i = 1:size(angcc_ov,1)                                                  %cortical index
    for j = 1:size(angcc_ov,2)                                              %cortical index
        if i ~= j                                                           %do not compare a region with itself
            disp(i); disp(j);

            % Ipsilateral or Contralateral
            cort_hemi_i = cpatch(i).hemi;                                   %ith patch hemisphere
            cort_hemi_j = cpatch(j).hemi;                                   %jth patch hemisphere
            if strcmp(cort_hemi_i, cort_hemi_j)
                same_hemisph(i,j) = 1;                                      %both cort in same hemisphere
                all_ang_same = [all_ang_same angcc_ov{i,j}(:)'];            %all the angles for ipsilateral
            elseif ~strcmp(cort_hemi_i, cort_hemi_j)
                same_hemisph(i,j) = 2;                                      %both cort in opposite hemispheres
                all_ang_opp = [all_ang_opp angcc_ov{i,j}(:)'];              %all the angles for contralateral
            end
    
            % Angle Distributions
            meanov_ang(i,j) = mean(angcc_ov{i,j});                          %mean principal angle for all representative modes
            if ~isempty(angcc_ov{i,j})
                minov_ang(i,j) = min(angcc_ov{i,j});                        %min principal angle for all representative modes
                maxov_ang(i,j) = max(angcc_ov{i,j});                        %max principal angle for all representative modes
            else
                minov_ang(i,j) = nan;                                       %empty
                maxov_ang(i,j) = nan;                                       %empty
            end
            %num_lowang(i,j) = length(find(angcc_ov{i,j}<20)); %# angles that are lower than 20 degrees
        else
            same_hemisph(i,j) = 0;                                          %indicates same region
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
            aopp = minov_ang(same_hemisph==2);      aopp = aopp(~isnan(aopp));
            savename = strcat('cc_',chtype,'_min_principal_angles');
        case 1
            asame = meanov_ang(same_hemisph==1);    asame = asame(~isnan(asame));    
            aopp = meanov_ang(same_hemisph==2);     aopp = aopp(~isnan(aopp));
            savename = strcat('cc_',chtype,'_mean_principal_angles');
        case 2
            asame = maxov_ang(same_hemisph==1);     asame = asame(~isnan(asame));
            aopp = maxov_ang(same_hemisph==2);      aopp = aopp(~isnan(aopp));
            savename = strcat('cc_',chtype,'_max_principal_angles');
        case 3
            asame = all_ang_same;       
            aopp = all_ang_opp;
            savename = strcat('cc_',chtype,'_all_principal_angles');
        case 4
            a_all = [all_ang_same all_ang_opp];
            savename = strcat('cc_',chtype,'_all_principal_angles_nohemis');
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