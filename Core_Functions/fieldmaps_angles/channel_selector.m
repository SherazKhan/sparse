function [chnames, ch_sel] = channel_selector(meastype, meas, chtype)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Select channels corresponding to meastype in an evoked measurement 
% Evoked measurement should be read in with MNE MATLAB Toolbox 
% Written by Pavitra Krishnaswamy
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    include = {}; exclude = meas.info.bads;

    if strcmp(meastype,'meg')                               
    
        want_meg = true; want_eeg = false; want_stim = false; 
        ch_sel = fiff_pick_types(meas.info,want_meg,want_eeg,want_stim,include,exclude);
        chnames = {meas.info.ch_names{ch_sel}};
        coil_type = [meas.info.chs(:).coil_type];           % coils differentiate mags vs grads   
        coil_type = coil_type(1:length(ch_sel));            % concatenate to length of #channels
        cij = ismember(coil_type, [3022 3024]);             % mag indexing by cij
        if strcmp(chtype, 'mag')
            ch_sel = ch_sel(cij);       chnames = chnames(cij);
        elseif strcmp(chtype, 'grad')
            ch_sel = ch_sel(~cij);      chnames = chnames(~cij);
        else
            disp('all channels');
            %ch_sel = ch_sel;            chnames = chnames;
        end
        
    elseif strcmp(meastype, 'eeg') && strcmp(chtype, 'eeg')
        
        want_meg = false; want_eeg = true; want_stim = false;
        ch_sel = fiff_pick_types(meas.info,want_meg,want_eeg,want_stim,include,exclude);
        chnames = {meas.info.ch_names{ch_sel}};   
    
    end
    
end