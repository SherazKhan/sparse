function plot_fcc_allpairs(angfcc_ov, cpatch, saveon, chtype, choice_fn)

%% Collate Principal Angle Distributions
all_ang = [];
for i = 1:length(cpatch)
    minov_ang(i) = min(angfcc_ov{i});                                       %min principal angle for all representative modes
    meanov_ang(i) = mean(angfcc_ov{i});                                     %mean principal angle for all representative modes
    maxov_ang(i) = max(angfcc_ov{i});                                       %max principal angle for all representative modes
    all_ang = [all_ang; angfcc_ov{i}];                                      %all the angles
end

%% Plot Histograms
for choose_fn = choice_fn:1:choice_fn
    % Collate Min, Mean, Max or All Angles
    switch choose_fn 
        case 0
            choose_ang = minov_ang;     savename = strcat('fcc_',chtype,'_min_principal_angles');
        case 1
            choose_ang = meanov_ang(~isnan(meanov_ang));    
            savename = strcat('fcc_',chtype,'_mean_principal_angles');
        case 2
            choose_ang = maxov_ang;     savename = strcat('fcc_',chtype,'_max_principal_angles');
        case 3
            choose_ang = all_ang;       savename = strcat('fcc_',chtype,'_all_principal_angles');
    end

    % Plot Histograms of Collated Angles
    figure, set(gcf,'color','white');
    %X = 0:2:90;N = hist(choose_ang, X); bar(X, N./sum(N)); xlim([0 90]);
    [~,X] = hist(choose_ang); histnorm(choose_ang,X); %hist(choose_ang); 
    xlabel('Principal Angle (Degrees)','FontSize',14); ylabel('Normalized Histogram','FontSize',14); 
    title('Principal Angles','FontSize',14); set(gca,'FontSize',14);    
    clear choose_ang X;
    
    % Save Histograms of Collated Angles
    if saveon
        print(gcf,'-depsc',savename);
    end
end

end