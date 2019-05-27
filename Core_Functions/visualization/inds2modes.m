function [hipsurf_ind1, reg_index, curr_str, curr_str_greedy, cpatch, sdiv_fwd_plot, all_regions] = inds2modes(sdiv_fwd, ...
         regname, all_regions, lk, trimP, ind_sfwd, ssubdiv, svolume, currstr, snummodes, clkmax, cp)

%% Index all Variables by # Modes
sdiv_srcspace0 = sum(cell2mat(ind_sfwd'));                                  % sdiv_fwd indices in joint source space
sdiv_fwd_plot0 = sdiv_fwd(sdiv_srcspace0 == 1);                             % sdiv_fwd of interest to us
sdiv_fwd_plot = []; sdiv_srcspace = [];                                     % initialize versions to account # modes
for i = 1:length(sdiv_fwd_plot0)
    if snummodes == 1
      sdiv_fwd_plot = [sdiv_fwd_plot; sdiv_fwd_plot0(i)];
      sdiv_srcspace = [sdiv_srcspace; sdiv_srcspace0(i)];
    else
      sdiv_fwd_plot = [sdiv_fwd_plot; repmat(sdiv_fwd_plot0(i), sdiv_fwd_plot0(i).numwhite_modes, 1)];
      sdiv_srcspace = [sdiv_srcspace; repmat(sdiv_srcspace0(i), sdiv_fwd_plot0(i).numwhite_modes, 1)];
    end
end
sdiv_fwd_plot = sdiv_fwd_plot';                                             % to makeit consistent with size of sdiv_fwd_plot0
clear all_regions;
for k = 1:length(sdiv_fwd_plot)
    all_regions{k} = sdiv_fwd_plot(k).reg_name;                             % names of regions of interest to us
end
ind_hipsurfs = strfind(all_regions, 'hipsurf');                             % indicate when region names are hipsurf
hipsurf_ind1 = find(cellfun('isempty', ind_hipsurfs) == 0,1);               % first element of sdiv_fwd_plot that is hipsurf 

%% Indices for Cortical and Subcortical Regions - to Find # Modes
cnt = 1;
for i = 1:(sum(lk<=clkmax))/cp
	cpatch{i,1} = lk(cnt:cnt+cp-1); cnt = cnt+cp;
end
clear cnt;
count = 1;
for cc = 1:sum(lk<=clkmax)
    reg_index(count,:) = [lk(cc),1];
    count = count + 1;
end
for sreg = 1:length(regname)
	for sdiv = 1:length(ssubdiv{sreg})
        nummodes = length(svolume{sreg,sdiv});
    	reg_index(count:count+nummodes-1,:) = repmat([sreg,sdiv],nummodes,1);
        count = count + nummodes;
	end
end
clear count nummodes sreg sdiv;
curr_str = currstr;                                                         % all current strengths
curr_str_greedy = currstr(trimP);                                           % selected current strengths for storage
%figure, plot(curr_str); figure, plot(curr_str_greedy);                     % plot current strengths

end