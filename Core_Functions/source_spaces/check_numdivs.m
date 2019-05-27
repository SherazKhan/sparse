function [ssing_mean, s_currstr_mean, n, split_ind, split_srclocs, geom_currstr] = ...
    check_numdivs(srclocs, n_init, savefigname, ploton, rawG, c_currstr_mean, ...
    sreg_cd, volume, ndipoles, geom_currstr, perc_match, cov_currstr_lim, megW, meg_sel, varargin)
% Written by Pavitra Krishnaswamy

n = n_init; cont_flag = 1; count = 0; num_count = 15; cs_match = 0.2;	    %Set Parameters
while(cont_flag > 0)
    count = count + 1; disp(count);
    
    %% Divide by Specified Number of Divsiions
    disp(strcat('trying to divide by ... ', num2str(n)));
    [split_ind, split_srclocs] = split_subdiv(srclocs, n, savefigname, ploton);  
    
    %% Check Singular Values for this Case
    for j = 1:length(split_ind)
        split_rawG = rawG(:,split_ind{j});                                  %Split the Raw Forward Solution
        if nargin == 14
		split_whiteG = megW*split_rawG(meg_sel,:);                  %Whiten the Split Forward Matrix
        else
		eegW = varargin{1}; eeg_sel = varargin{2}; 
		split_whiteG = cat(1, megW*split_rawG(meg_sel,:), eegW*split_rawG(eeg_sel,:)); %Whiten the Split Forward Matrix and Concat
        end
        whiteS{j} = svd(split_whiteG,0);                                    %Compute SVD
        max_whiteS(j) = max(whiteS{j});                                     %Max Singular Value = Norm 
        mean_whiteS(j) = mean(whiteS{j});                                   %Mean Singular Value Across Modes
        realvol(j) = volume/ndipoles*floor(size(split_ind{j},2)/3);         %Physical volume of Split
    end

    %% Compute Representative Current Strengths and Compare with Cortical Current Strength
    ssing_mean = mean(max_whiteS);                                          %Representative Singular Value
    s_geom_currstr = sreg_cd*volume/n;
    s_currstr_mean = sreg_cd*volume/n*ssing_mean;                           %Representative Current Strength
    norm_diff_mean = s_currstr_mean/c_currstr_mean - 1;                     %Ratio with Cortical Current Strength
    disp(strcat('subcortical current strength minus cortical strength is ...', num2str(norm_diff_mean*100),'%'));
    
    %% Compute Variability in Current Strength Across Divisions
    currstr = sreg_cd*realvol;                                              %real current strengths of each subdivision
    cov_currstr = std(currstr)/mean(currstr);                               %coefficient of variation in current strength
    disp(strcat('subcortical current strength variability is ...', num2str(cov_currstr*100),'%'));    
    
    %% Compare with Other Subcortical Geometric Current Strengths
    disp(strcat('geometric strength: ...', num2str(s_geom_currstr)));       %geometric current strength at present
    diff_geom_currstr = s_geom_currstr/mean(geom_currstr) - 1               %how different is the current strength from other regions
    
    %% Update Number of Divisions and Decide on Iteration
    if abs(norm_diff_mean) > perc_match  || abs(diff_geom_currstr) > cs_match%|| cov_currstr > 0.15 %|| abs(norm_diff_max) > 0.2
        cont_flag = 1;                                                      %Must Continue Iteration
        if norm_diff_mean < 0                                               %Dividing too much
            n = n - 2;
        elseif norm_diff_mean > 0 || abs(diff_geom_currstr) > cs_match 	    %Dividing too little (can also || cov_currstr > 0.15)
            n = n + 2;
        end
    else
        cont_flag = 0;                                                      %Stop the Divisions
        n = length(split_ind);                                              %Actual # Divisions
    end
    
    %% Not Converging: Relax Criteria a Bit
    if count > num_count && cont_flag == 1
        if num_count < 25
            perc_match = perc_match*2                                       %Not Converging: Relax the Cortical = Subcortical Constraint
            num_count = num_count + 10;                                     %Up the num count to give it a significant chance to converge again
        else
            disp('over 25');
            if cov_currstr > cov_currstr_lim
            	n = round(sreg_cd*volume/mean(geom_currstr));               %set number of divisions and be done, it will not converge
                perc_match = 0.5;       
            else
                perc_match = perc_match*2;
                cs_match = cs_match*2;
            end
            num_count = num_count + 10; 
          
        end
    end
end

if strcmp(savefigname(end-1:end),'lh')
    geom_currstr = s_geom_currstr;
else
    geom_currstr = [geom_currstr s_geom_currstr];
end
    
end