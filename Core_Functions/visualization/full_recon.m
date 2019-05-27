function [Xc_plot] = full_recon(alg,ico,fwd,p,patchno,prefix,source,stc,subdiv,trimP,X_hat,L,curr_den_fname,scalarX, varargin)
%function full_recon(alg,ico,fwd,p,patchno,prefix,source,stc,subdiv,trimP,X_hat,L, varargin)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% RECONSTRUCTION OF SOURCE ESTIMATES OVER FULL-ORDER SOURCE SPACE 
% Written by Gabriel Obregon-Henao (obregon@nmr.mgh.harvard.edu)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Initialize
T = size(X_hat,2);
Xc_plot(1:length(patchno)) = struct('numdipoles',[],'cstrdip',[],'mode',[],'orient',[],'dip_or',[],'dip',[],...
        'dip_res',[],'dip_resor',[],'dip_resovern',[],'dip_resnorm',[]);
Xfdip = zeros(fwd.nsource,T);
Xfdip_resovern = zeros(fwd.nsource,T);
load(curr_den_fname,'cort_cd');                                             %cortical current strength

%% Project Estimates to Full-order Source Space
for ii = 1:length(patchno)
    % Basic Projection of Amplitudes to Dipole Space
    V = ico(subdiv).patch(patchno(ii)).V(:,1:p);                            %area of chosen spatial modes
    A = ico(subdiv).patch(patchno(ii)).A;                                   %area of chosen patch
    Xc_plot(ii).numdipoles = size(V,1);                                     %# dipoles
    Xc_plot(ii).cstrdip = A*cort_cd/Xc_plot(ii).numdipoles;                 %current strength per dipole
    Xc_plot(ii).mode = X_hat(p*trimP(ii)-(p-1):p*trimP(ii),:);              %mode amplitudes - all modes of selected patch
    Xfdip(ico(subdiv).patch(patchno(ii)).lk,:) = V*Xc_plot(ii).mode;        %dipole amplitudes
    Xc_plot(ii).dip = Xfdip(ico(subdiv).patch(patchno(ii)).lk,:);           %dipole amplitudes
    
    % Compute Resultant Dipole Time Series Estimates in Hypersampled Source Space
    Xc_plot(ii).orient = fwd.source_nn(ico(subdiv).patch(patchno(ii)).lk,:);%orientations of dipoles - MEG coordframe
    nn(:,:,1) = Xc_plot(ii).orient; nn = repmat(nn,[1,1,T]);
    Xc_plot(ii).dip_or(:,:,1) = squeeze(nn(:,1,:)).*Xc_plot(ii).dip;        %1 X #timepoints
    Xc_plot(ii).dip_or(:,:,2) = squeeze(nn(:,2,:)).*Xc_plot(ii).dip;        %1 X #timepoints
    Xc_plot(ii).dip_or(:,:,3) = squeeze(nn(:,3,:)).*Xc_plot(ii).dip;        %1 X #timepoints
    [Xc_plot(ii).dip_res, Xc_plot(ii).dip_resor] = vec_sum(Xc_plot(ii).dip_or);   %net activity over all dipoles
    Xc_plot(ii).dip_resovern = Xc_plot(ii).dip_res/Xc_plot(ii).numdipoles;  %average activity per dipole
    Xc_plot(ii).dip_resnorm = Xc_plot(ii).dip_res/(Xc_plot(ii).cstrdip*Xc_plot(ii).numdipoles); %average activity per dipole (unitless) 
    Xfdip_res(ico(subdiv).patch(patchno(ii)).lk,:) = repmat(Xc_plot(ii).dip_res,Xc_plot(ii).numdipoles,1); %average amplitude/dipole
    Xfdip_resovern(ico(subdiv).patch(patchno(ii)).lk,:) = repmat(Xc_plot(ii).dip_resovern,Xc_plot(ii).numdipoles,1); %average amplitude/dipole
    clear nn 
end

%% Compute Number of Patches in Left Hemisphere
puse = length(source(subdiv).hemisph(1).pinfo);
p1 = patchno <= puse;

%% Create STC File Structure Fields for Left Hemisphere - Amplitude of Resultant Across Dipoles
if sum(p1) > 0
    disp('left');
    if nargin == 14
      lk = [ico(subdiv).patch(patchno(p1)).lk];
    else
      patchnum = varargin{1};
      lk = [ico(subdiv).patch(patchnum(p1)).lk];
    end
    stc(1).vertices = fwd.src(1).vertno(lk);
    stc(1).data = Xfdip_res(lk,:);                                          %scalarX if Xfdip_resovern
else 
    stc(1).vertices = [];   stc(1).data = [];                               %Initialize for Results View
end

%% Create STC File Structure Fields for Right Hemisphere - Amplitude of Resultant Across Dipoles
if sum(p1) < length(p1)            
    disp('right');
    if nargin == 14
      lk2 = [ico(subdiv).patch(patchno(~p1)).lk];
    else
      patchnum = varargin{1};
      lk = [ico(subdiv).patch(patchnum(~p1)).lk];
    end
    stc(2).vertices = fwd.src(2).vertno(lk2-length(fwd.src(1).vertno));
    stc(2).data = Xfdip_res(lk2,:);                                         %scalarX if Xfdip_resovern
else 
    stc(2).vertices = [];   stc(2).data = [];                               %Initialize for Results View
end

%% Write STC Files
switch alg
    case 'MNE'
        mne_write_stc_file1(strcat(prefix,'_ico_',num2str(subdiv),'p_SPIGH_L1_',num2str(L),'_patchno_',num2str(length(patchno)),'_est-lh.stc'),stc(1))
        mne_write_stc_file1(strcat(prefix,'_ico_',num2str(subdiv),'p_SPIGH_L1_',num2str(L),'_patchno_',num2str(length(patchno)),'_est-rh.stc'),stc(2))
    case 'dSPM'
        mne_write_stc_file1(strcat(prefix,'_ico_',num2str(subdiv),'p_SPIGH-dSPM_L1_',num2str(L),'_patchno_',num2str(length(patchno)),'_est-lh.stc'),stc(1))
        mne_write_stc_file1(strcat(prefix,'_ico_',num2str(subdiv),'p_SPIGH-dSPM_L1_',num2str(L),'_patchno_',num2str(length(patchno)),'_est-rh.stc'),stc(2)) 
    case 'sMAP-EM'
        mne_write_stc_file1(strcat(prefix,'_ico_',num2str(subdiv),'p_SPIGH-EM_L1_',num2str(L),'_patchno_',num2str(length(patchno)),'_est-lh.stc'),stc(1))
        mne_write_stc_file1(strcat(prefix,'_ico_',num2str(subdiv),'p_SPIGH-EM_L1_',num2str(L),'_patchno_',num2str(length(patchno)),'_est-rh.stc'),stc(2))        
end

end