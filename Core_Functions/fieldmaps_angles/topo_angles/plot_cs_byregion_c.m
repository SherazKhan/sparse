function clikely_angle = plot_cs_byregion_c(angcs_ov, sdiv, cpatch, saveon, ang_flag)

%% Initialize
call_ang_same = cell(1,size(angcs_ov,2));       
call_ang_opp = cell(1,size(angcs_ov,2));        
clikely_angle = zeros(1,size(angcs_ov,2));      

for j = 1:size(angcs_ov,2)                                                  %cortical index
    %% Collate Principal Angle Distributions
    for i = 1:size(angcs_ov,1)                                              %subcortical index
        disp(i); disp(j);
        
        % Ipsilateral or Contralateral
        subc_hemi = sdiv(i).hemi;                                           %subcortical hemisphere
        cort_hemi = cpatch(j).hemi;                                         %cortical hemisphere
        if strcmp(subc_hemi, cort_hemi)
            same_hemisph(i,j) = 1;                                          %cort, scort in same hemisphere
            call_ang_same{j} = [call_ang_same{j} angcs_ov{i,j}(:)'];        %all the angles for ipsilateral
        elseif ~strcmp(subc_hemi, cort_hemi)
            same_hemisph(i,j) = 0;                                          %cort, scort in opposite hemispheres
            call_ang_opp{j} = [call_ang_opp{j} angcs_ov{i,j}(:)'];          %all the angles for contralateral
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
    
    %% Plot Histograms
    for choose_fn = 3:3                                                     % Min, Mean, Max or All Angles
        %pause(1); 
        close; 
        figure(5+choose_fn), set(gcf,'color','white');
        switch choose_fn 
            case 0
                a = minov_ang(:,j);
                savename = strcat('min');
            case 1
                a = meanov_ang(:,j);
                savename = strcat('mean');
            case 2
                a = maxov_ang(:,j);    
                savename = strcat('max');
            case 3
                a = sort([call_ang_same{j} call_ang_opp{j}]);
                savename = strcat('all');
        end

        if strcmp(ang_flag, 'mode')
          [ns,Xs] = hist(a); histnorm(a,Xs);
          %[ns,Xs] = hist(a,20); histnorm(a,Xs); 
          [~, maxns] = max(ns);   clikely_angle(j) = Xs(maxns);
          if isnan(clikely_angle(j))
            clikely_angle(j) = 0;
          end
          title(strcat('Region',num2str(j),'Vs. All Other')); 
          xlabel('Principal Angle (Degrees)','FontSize',14); ylabel('Normalized Histogram','FontSize',14); 
          xlim([0 90]); set(gca,'FontSize',14);    
          clear a Xs;
          if saveon
            saveas(gcf,strcat(pwd,'/cort_vs_scort_hists/cpatch',num2str(j),'_vs_alls_',savename,'ang.jpg'));
          end
        elseif strcmp(ang_flag,'median')
          clikely_angle(j) = median(a);
        end
    end
end

end