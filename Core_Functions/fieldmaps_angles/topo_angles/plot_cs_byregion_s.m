function slikely_angle = plot_cs_byregion_s(angcs_ov, sdiv, cpatch, saveon, ang_flag)

%% Initialize
sall_ang_same = cell(1,size(angcs_ov,1));       
sall_ang_opp = cell(1,size(angcs_ov,1));        
slikely_angle = zeros(1,size(angcs_ov,1));      

for i = 1:size(angcs_ov,1)                                                  %subcortical index
    %% Collate Principal Angle Distributions
    for j = 1:size(angcs_ov,2)                                              %cortical index
        disp(i); disp(j);
        
        % Ipsilateral or Contralateral
        subc_hemi = sdiv(i).hemi;                                           %subcortical hemisphere
        cort_hemi = cpatch(j).hemi;                                         %cortical hemisphere
        if strcmp(subc_hemi, cort_hemi)
            same_hemisph(i,j) = 1;                                          %cort, scort in same hemisphere
            sall_ang_same{i} = [sall_ang_same{i} angcs_ov{i,j}(:)'];        %all the angles for ipsilateral
        elseif ~strcmp(subc_hemi, cort_hemi)
            same_hemisph(i,j) = 0;                                          %cort, scort in opposite hemispheres
            sall_ang_opp{i} = [sall_ang_opp{i} angcs_ov{i,j}(:)'];          %all the angles for contralateral
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
                a = minov_ang(i,:);
                savename = strcat('min');
            case 1
                a = meanov_ang(i,:);
                savename = strcat('mean');
            case 2
                a = maxov_ang(i,:);    
                savename = strcat('max');
            case 3
                a = sort([sall_ang_same{i} sall_ang_opp{i}]);
                savename = strcat('all');
        end
        
        if strcmp(ang_flag, 'mode')
          [ns,Xs] = hist(a); histnorm(a,Xs);
          %[ns,Xs] = hist(a,20); histnorm(a,Xs); 
          [~, maxns] = max(ns);   slikely_angle(i) = Xs(maxns);
          if isnan(slikely_angle(i))
            slikely_angle(i) = 0;
          end
          title(strcat('Region',num2str(i),'Vs. All Other')); 
          xlabel('Principal Angle (Degrees)','FontSize',14); ylabel('Normalized Histogram','FontSize',14); 
          xlim([0 90]); set(gca,'FontSize',14);    
          clear a Xs;
          if saveon
            saveas(gcf,strcat(pwd,'/scort_vs_cort_hists/sdiv',num2str(i),'_vs_allc_',savename,'ang.jpg'));
          end
        elseif strcmp(ang_flag,'median')
          slikely_angle(i) = median(a);
        end
    end   
end

end