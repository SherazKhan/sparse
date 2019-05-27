function [Xls_plot_cond, varargout] = roi_sum(rois, unique_rois, uni_divs, Xls_plot_alg, alg_display, varargin)

disp(alg_display);
ll = length(unique_rois);
Xls_plot(1:ll) = struct('equiv_dips',[]);
Xls_plot_cond(1:ll) = struct('numdipoles',[],'cstrdip',[],'orient',[],'dip_or',[],'dip_res',[],'dip_resor',[],'dip_resovern',[],'dip_resnorm',[],'mode_indices',[],'mode',[]);
if nargin > 5
    rois_revise = varargin{1};
end

for i = 1:length(unique_rois)
    roi_patches = uni_divs(strcmp(rois, unique_rois(i)))
    if length(roi_patches) > 1
        for j = 1:length(roi_patches)
            Xls_resor(j,:,:) = Xls_plot_alg(roi_patches(j)).dip_resor';         %nX3XT
            Xls_resmag(j,:) = Xls_plot_alg(roi_patches(j)).dip_res;             %nXT
        end
        Xls_plot(i).equiv_dips(:,:,1) = squeeze(Xls_resor(:,1,:)).*Xls_resmag;  %orientation 1
        Xls_plot(i).equiv_dips(:,:,2) = squeeze(Xls_resor(:,2,:)).*Xls_resmag;  %orientation 2
        Xls_plot(i).equiv_dips(:,:,3) = squeeze(Xls_resor(:,3,:)).*Xls_resmag;  %orientation 3
        [Xls_plot_cond(i).dip_res, Xls_plot_cond(i).dip_resor] = vec_sum(Xls_plot(i).equiv_dips); 
        if nargin > 5
            rois_revise(roi_patches) = strcat(unique_rois(i),'-all'); 
        end
    else
        Xls_plot_cond(i).dip_res = Xls_plot_alg(roi_patches).dip_res;
        Xls_plot_cond(i).dip_resor = Xls_plot_alg(roi_patches).dip_resor;
    end
end
if nargin > 5
    varargout{1} = rois_revise;
end
      
end