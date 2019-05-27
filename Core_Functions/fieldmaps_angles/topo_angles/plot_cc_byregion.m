function likely_angle = plot_cc_byregion(angcc_ov, cpatch, saveon, ang_flag)

%% Initialize
all_ang_same = cell(1,size(angcc_ov,1)); 
all_ang_opp = cell(1,size(angcc_ov,1));
likely_angle = zeros(1,size(angcc_ov,1));

for i = 1:size(angcc_ov,1)                                                  %cortical index
    %% Collate Principal Angle Distributions
    for j = 1:size(angcc_ov,2)                                              %cortical index
        if i ~= j                                                           %do not compare a region with itself
            disp(i); disp(j);

            % Ipsilateral or Contralateral
            cort_hemi_i = cpatch(i).hemi;                                   %ith patch hemisphere
            cort_hemi_j = cpatch(j).hemi;                                   %jth patchhemisphere
            if strcmp(cort_hemi_i, cort_hemi_j)
                same_hemisph(i,j) = 1;                                      %both cort in same hemisphere
                all_ang_same{i} = [all_ang_same{i} angcc_ov{i,j}(:)'];      %all the angles for ipsilateral
            elseif ~strcmp(cort_hemi_i, cort_hemi_j)
                same_hemisph(i,j) = 2;                                      %both cort in opposite hemispheres
                all_ang_opp{i} = [all_ang_opp{i} angcc_ov{i,j}(:)'];        %all the angles for contralateral
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
            saveas(gcf,strcat(pwd,'/cort_vs_cort_hists/cpatch',num2str(i),'_vs_allc_',savename,'ang.jpg'));
          end
        elseif strcmp(ang_flag,'median')
          likely_angle(i) = median(a);
        end
    end
end

end