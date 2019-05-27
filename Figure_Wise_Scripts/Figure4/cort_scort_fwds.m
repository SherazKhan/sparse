function [Gtheta, Gn, currstr, lk_flag, lk, svolume, clkmax, reg_Gtheta] = cort_scort_fwds(lk, ico, cp, cort_cd, ...
sdiv_fwd, sregname, sregnums, nummodes, choose_scort)

if ~isempty(lk)
    %% Scale and Merge Whitened Cortical Forward Solutions
    Gcn = []; Gctheta = []; clk = []; %A = []; 
    if nummodes
        cp = 1;                                                             %choose only 1 top mode per patch
    end
    for ii = 1:length(lk)                                                   %for each cortical patch
        cpatch(ii) = ico.patch(lk(ii));                                     %read in cortical patch
        %if cpatch(ii).exclude == 0                                          %don't exclude as patch not on medial wall
        Gcn = cat(2,Gcn,cpatch(ii).U(:,1:cp));                              %modes
        Gctheta = cat(2,Gctheta,cpatch(ii).U(:,1:cp)*cpatch(ii).S(1:cp,1:cp)); %forward soln
        %A = [A repmat(cpatch(ii).A,1,cp)];                                 %individual patch areas
        clk = [clk repmat(lk(ii),1,cp)];                                    %as each mode is independent, repeat patchno for each mode
        %end
    end
    A = repmat(mean([ico.patch(:).A]),1,length(clk));                       %average area of any cortical patch
    c_curr_str = A*cort_cd;                                                 %average current strength for any cortical patch
    Gn = Gcn;                                                               %U's dont scale by current strength
    Gtheta = Gctheta*diag(c_curr_str);                                      %scale Gctheta by cortical current strength
    currstr = c_curr_str;                                                   %current strength
    lk = clk;                                                               %cortical patches to go into subcortical soln
    lk_flag = ones(1,length(lk));                                           %flag to cortical indicators (1)
    clkmax = length(ico.patch);%max(max(lk), length(lk_flag));              %maximum of all cortical patches while accounting for modes

else
    %% Empty Cortical Forward Solutions
    Gn = []; Gtheta = []; currstr = []; lk = []; lk_flag = []; clkmax = 0; 
    disp('no cortical forward solutions');
end

%% Read in Subcortical Fwd Structure For Easy Indexing
temp_lkmax = clkmax + 1; lkmax = clkmax;                                    %initialize
for k = 1:length(sdiv_fwd)                          
    if nummodes
    sp(k) = 1;                                                              %#white modes 
    else
    sp(k) = sdiv_fwd(k).numwhite_modes;                                     %#white modes
    end
    all_regions{k} = sdiv_fwd(k).reg_name;                                  %names of regions
    sfwd(k).currstr = repmat(sdiv_fwd(k).currstr,1,sp(k));                  %current strength
    sfwd(k).Gstheta_p = sdiv_fwd(k).whiteGstheta(:,1:sp(k));                %whitened G cols (scaled by currstr)
    sfwd(k).Us_p = sdiv_fwd(k).whiteUs(:,1:sp(k));                          %whitened modes (not scaled by currstr)
    sfwd(k).index = repmat(temp_lkmax, 1, sp(k));                           %repeat region index
    temp_lkmax = sfwd(k).index(1) + 1;                                      %update starting point for region indices
    %temp_lkmax = temp_lkmax + length(sfwd(k).index) - 1;                   %update starting point for region indices
end
left_right = lk < length(ico.patch)/2;                                      %if left this is 1
reg_Gtheta = cell(1, length(lk)); 
left = repmat({'lcortical'},1,length(lk));    right = repmat({'rcortical'},1,length(lk));
reg_Gtheta(left_right) = left(left_right);    reg_Gtheta(~left_right) = right(~left_right); %initialize region names

if ~strcmp(choose_scort,'none')
    %% Scale and Merge Whitened Subcortical Forward Solutions
    for i = 1:length(sregnums)
        %Identify ith Region
        regname{i} = sregname{sregnums(i)};                                 %name of ith region selected
        ind_sfwd = strncmp(regname{i}, all_regions, 4);                     %find indices of sdiv_fwd that match regname{i}
        ind = find(ind_sfwd); 
        
        %Load in Fwd Solutions of ith Region
        Gstheta = [sfwd(ind_sfwd).Gstheta_p];                               %all the Gstheta for ith region
        Gsn = [sfwd(ind_sfwd).Us_p];                                        %all the Us for ith region
        Gtheta = cat(2,Gtheta, Gstheta);                                    %put Gstheta into overall Gtheta
        Gn = cat(2, Gn, Gsn);                                               %put Gsn into overall Gn
        currstr = [currstr,sfwd(ind_sfwd).currstr];                         %put sfwd.currstr into overall currstr
        reg_Gtheta = [reg_Gtheta repmat(all_regions(ind(1)), 1, size(Gstheta,2))];
        
        %Define Division Indices for ith Region
        ssubdiv{i} = find(ind_sfwd);                                        %subdivisions of region i
        snummodes = sp(ind_sfwd);                                           %#modes to use in region i's subdivisions
        c = 1;                                                              %counter
        for j = 1:length(ssubdiv{i})                                        %length(ssubdiv{i}) is #subdivisions
            spj = snummodes(j);                                             %#modes for ith region jth subdivision
            slk = repmat(lkmax+1, 1, spj);                                  %lk indices for ith region
            svolume{i,j} = slk;                                             %indices for ith region jth subdivision
            lk = [lk slk]; lkmax = max(lk);                                 %Add to Global Pool, Index Accordingly
            c = c + spj;
        end
        %svolume{i} = slk;                                                  %agnostic to subdiv index 
        lk_flag = [lk_flag zeros(1,size(Gsn,2))];                           %flag to subcortical indicators (0)
    end
else
    %% Empty Subcortical Forward Solutions
    svolume = [];
    disp('no subcortical forward solutions'); 
end

end