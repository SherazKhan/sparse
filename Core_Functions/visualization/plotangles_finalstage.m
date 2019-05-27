function plotangles_finalstage(lk, clkmax, Gtheta, c_curr_str, hipsurf_ind1, all_regions)

% Indices of Various Parts of Gtheta
num_cort = length(c_curr_str);
num_lhip = sum(strcmp(all_regions,'lhipsurf'));
num_rhip = sum(strcmp(all_regions,'rhipsurf'));
cortind = 1:num_cort;                                                       %cortical
firsthip = hipsurf_ind1+num_cort;  ind_left_hip = firsthip:firsthip+num_lhip-1; %left hippocampal
secondhip = firsthip+num_lhip; ind_right_hip = secondhip: secondhip+num_rhip-1; %right hippocampal 
ind_bs = num_cort + find(strcmp(all_regions,'bsred'));                      %brainstem
left_right_trind = num_cort + find(strcmp(all_regions,'ra'),1,'first');     %right volumes
left_vols_ind = (num_cort +1):(left_right_trind-1);                         %left volumes
right_vols_ind = + [left_right_trind:min([ind_bs(1), ind_left_hip(1), ind_right_hip(1)])-1];

% Rearrange columns of Gtheta for Easy Reading of Angle Matrices
newGtheta = Gtheta(:,cortind); %size(newGtheta)                             %cortical surface
newGtheta = cat(2,newGtheta, Gtheta(:,ind_left_hip)); %size(newGtheta)      %left hippocampal surface
newGtheta = cat(2,newGtheta, Gtheta(:,left_vols_ind)); %size(newGtheta)     %left volumes
newGtheta = cat(2, newGtheta, Gtheta(:,ind_bs)); %size(newGtheta)           %brainstem volumes
newGtheta = cat(2, newGtheta, Gtheta(:,ind_right_hip)); %size(newGtheta)    %right hippocampal surface
newGtheta = cat(2, newGtheta, Gtheta(:,right_vols_ind)); %size(newGtheta)   %right volumes

% Cortical vs. Subcortical Mode Pairs
x = sum(lk<=clkmax);
for i = 1:x
    for j = x+1:size(newGtheta,2)
        angcs(i,j-x) = subspacea(newGtheta(:,i), newGtheta(:,j))*180/pi;
    end
end
figure, set(gcf,'color','white'); imagesc(angcs); axis xy; set(gca,'FontSize',14);
caxis([0 90]); c = colorbar; ylabel(c,'Angles (Degrees)','FontSize',14);
ylabel('Cortical Modes','FontSize',14); xlabel('Subcortical Modes','FontSize',14); 
title('Angles between Cortical and Subcortical Modes in Final Stage','FontSize',14);

% Subcortical vs. Subcortical Mode Pairs
for i = x+1:size(newGtheta,2)
    for j = x+1:size(newGtheta,2)
        angss(i-x,j-x) = subspacea(newGtheta(:,i), newGtheta(:,j))*180/pi;
    end
end
figure, set(gcf,'color','white'); imagesc(angss); axis xy; set(gca,'FontSize',14);
caxis([0 90]); c = colorbar; ylabel(c,'Angles (Degrees)','FontSize',14);
ylabel('Subcortical Modes','FontSize',14); xlabel('Subcortical Modes','FontSize',14); 
title('Angles between Subcortical Modes in Final Stage','FontSize',14);

end