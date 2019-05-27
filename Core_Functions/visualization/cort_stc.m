function cort_stc(prefix, meastype, datatype, lk, trimP, clkmax, cfwd, Xls_plot, stc, stc_name)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Generate STC movies of cortical sources estimated in joint deep + superficial pursuit hierarchy
% Written by Pavitra Krishnaswamy (pavitrak@nmr.mgh.harvard.edu)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Initialization
load(strcat(prefix,'_',meastype,'_plr_prewhite_',datatype,'_cSVDs.mat'), 'source','ico'); 
ico = ico(3); puse = length(source(3).hemisph(1).pinfo); clear source;
sel_lk = lk(trimP);

%% Create STC File Structure Fields
for i = 1:length(sel_lk)
    if sel_lk(i) < clkmax
       if sel_lk(i) < puse
	   lk1 = ico.patch(sel_lk(i)).lk;
	   stc(1).vertices = cfwd.src(1).vertno(lk1);
	   stc(1).data = repmat(Xls_plot(trimP(i)).dip_res, Xls_plot(trimP(i)).numdipoles, 1);
       else
	   lk2 = ico.patch(sel_lk(i)).lk - length(cfwd.src(1).vertno);
	   stc(2).vertices = cfwd.src(2).vertno(lk2);
	   stc(2).data = Xls_plot(trimP(i)).dip_res;
       end
    end
end

%% Write STC Files
mne_write_stc_file1(strcat(prefix,'_', stc_name,'-lh.stc'), stc(1));
mne_write_stc_file1(strcat(prefix,'_', stc_name,'-rh.stc'), stc(2));

end