function likely_angle = plot_ss_byregion(angss_ov, sdiv, saveon, ang_flag)

%% Initialize
all_ang_same = cell(1,size(angss_ov,1)); 
all_ang_opp = cell(1,size(angss_ov,1));
likely_angle = zeros(1,size(angss_ov,1));

for i = 1:size(angss_ov,1)                                                  %subcortical index
    %% Collate Principal Angle Distributions for ith Region 
    for j = 1:size(angss_ov,2)                                              %subcortical index
        if i ~= j                                                           %do not compare a region with itself
            disp(i); disp(j);

            % Ipsilateral or Contralateral
            subc_hemi_i = sdiv(i).hemi;                                     %ith subcortical hemisphere
            subc_hemi_j = sdiv(j).hemi;                                     %jth subcortical hemisphere
            if strcmp(subc_hemi_i, subc_hemi_j)
                same_hemisph(i,j) = 1;                                      %both scort in same hemisphere
                all_ang_same{i} = [all_ang_same{i} angss_ov{i,j}(:)'];      %all the angles for ipsilateral
            elseif ~strcmp(subc_hemi_i, subc_hemi_j)
                same_hemisph(i,j) = 2;                                      %both scort in opposite hemispheres
                all_ang_opp{i} = [all_ang_opp{i} angss_ov{i,j}(:)'];        %all the angles for contralateral
            end

            % Angle Distributions
            meanov_ang(i,j) = mean(angss_ov{i,j});                          %mean principal angle for all representative modes
            if ~isempty(angss_ov{i,j})
                minov_ang(i,j) = min(angss_ov{i,j});                        %min principal angle for all representative modes
                maxov_ang(i,j) = max(angss_ov{i,j});                        %max principal angle for all representative modes
            else
                minov_ang(i,j) = nan;                                       %empty
                maxov_ang(i,j) = nan;                                       %empty
            end
            %num_lowang(i,j) = length(find(angss_ov{i,j}<20)); %# angles that are lower than 20 degrees
        else
            same_hemisph(i,j) = 0;                                          %indicates same region
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
                a = sort([all_ang_same{i} all_ang_opp{i}]);
                savename = strcat('all');
        end
        
        if strcmp(ang_flag, 'mode')
          [ns,Xs] = hist(a); histnorm(a,Xs);
          %[ns,Xs] = hist(a,20); histnorm(a,Xs); 
          [~, maxns] = max(ns);   likely_angle(i) = Xs(maxns);
          if isnan(likely_angle(i))
            likely_angle(i) = 0;
          end
          title(strcat('Region',num2str(i),'Vs. All Other')); 
          xlabel('Principal Angle (Degrees)','FontSize',14); ylabel('Normalized Histogram','FontSize',14); 
          xlim([0 90]); set(gca,'FontSize',14);    
          clear a Xs;
          if saveon
            saveas(gcf,strcat(pwd, '/scort_vs_scort_hists/sdiv',num2str(i),'_vs_alls_',savename,'ang.jpg'));
          end
        elseif strcmp(ang_flag,'median')
          likely_angle(i) = median(a);
        end
    end
    
end

end