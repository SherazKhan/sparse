function [xls_resmag,xls_resor] = vec_sum(xls_dip_or)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%xls_dip_or should be #dipoles X #timepoints X #orientations
%xls_res is #1 X time points 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
%% Magnitude of the Sum
xls1 = sum(xls_dip_or(:,:,1),1);                                            %vector sum along 1st orientation
xls2 = sum(xls_dip_or(:,:,2),1);                                            %vector sum along 2nd orientation
xls3 = sum(xls_dip_or(:,:,3),1);                                            %vector sum along 3rd orientation
xls_resmag = (sqrt(xls1.^2 + xls2.^2 + xls3.^2));                           %net activity across all orientations
    
%% Orientation of the Dipole Sum
xls_resor(:,1) = xls1./xls_resmag;                                          %orientation vector as a function of time
xls_resor(:,2) = xls2./xls_resmag;                                          %orientation vector as a function of time
xls_resor(:,3) = xls3./xls_resmag;                                          %orientation vector as a function of time
    
end

