function [Xls_plot_modes, Xls_plot] = backproj_greedy(lk, clkmax, trimP, strimP, ctrimP, reg_index, svolume, cpatch, ...
         sdiv_fwd_plot, ico, Xls_scaled, currstr, hsrc, hfwd, csrc)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Project scaled estimates from eigenspace to physical space for plotting greedy estimates
% Written by Pavitra Krishnaswamy (pavitrak@nmr.mgh.harvard.edu)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Initialize Variables
Xls_plot(1:length(trimP)) = struct('numdipoles',[],'cstrdip',[],'mode',[],'orient',[],'dip_or',[],'dip',[],...
        'dip_res',[],'dip_resor',[],'dip_resovern',[],'dip_resnorm',[]);
T = size(Xls_scaled,2);                                                     % # time points
cnn = cat(1,csrc(1).nn,csrc(2).nn);                                         % all cortical orientations

for i = 1:length(trimP)
  s = strimP(i);                                                            % index wrt subcortical columns in joint srcspace
  c = ctrimP(i);                                                            % index wrt cortical columns in joint srscpace
  slkmax = sum(lk<=clkmax);                                                 % number of cortical columns in joint srcspace
  
  if s > 0                                                 
    %% Backproject Subcortical Scaled Estimates from Eigen Space to Physical Space
    sdiv = sdiv_fwd_plot(s-slkmax);                                         % sdiv_fwd relevant to trimP(i)
    mode_num = svolume{reg_index(s,1),reg_index(s,2)} == lk(s);             % index of chosen mode
    Vs = sdiv.whiteVs(:,mode_num);                                          % right eigenvector for chosen mode - num_dipoles X 1
    Xls_sdip{i} = Vs*Xls_scaled(s,:);                                       % backproject greedy modal estimates to dipoles
    Xls_plot(i).numdipoles = size(Vs,1)/3;                                  % number of dipoles
    Xls_plot(i).cstrdip = currstr(s)/Xls_plot(i).numdipoles;                % current strength per dipole
    Xls_plot(i).mode = Xls_scaled(s,:);                                     % mode current estimate
    if isempty(strfind(sdiv.reg_name,'hipsurf'))                            % subcortical volumes - so need to separate out 3 orients
        Xls_plot(i).orient = [1 1 1];
        Xls_plot(i).dip_or(:,:,1) = Xls_sdip{i}(1:3:end,:);                 % orientation 1
        Xls_plot(i).dip_or(:,:,2) = Xls_sdip{i}(2:3:end,:);                 % orientation 2
        Xls_plot(i).dip_or(:,:,3) = Xls_sdip{i}(3:3:end,:);                 % orientation 3
        Xls_plot(i).dip = sqrt(sum(Xls_plot(i).dip_or.^2,3));               % dipole current estimate (amplitude)
    else                                                                    % hippocampal surfaces - no need to separate out 3 orients
      Xls_plot(i).dip = Xls_sdip{i};                                        % dipole current estimate
      if isempty(strfind(sdiv.reg_name,'rh'));                              % left orientation vectors
          if isfield(sdiv, 'patchno')
          [~,~,hlk] = intersect([hsrc(1).pinfo{sdiv.patchno}],hfwd.src(1).vertno,'legacy'); 
          else
          [~,~,hlk] = intersect(hsrc(1).pinfo{sdiv.index},hfwd.src(1).vertno,'legacy'); 
          end
          Xls_plot(i).orient = hsrc(1).nn(hlk,:);
          %rr{i} = hsrc(1).rr(hlk,:);
      else                                                                  % right orientation vectors
          if isfield(sdiv, 'patchno')
          [~,~,hrk] = intersect([hsrc(2).pinfo{sdiv.patchno}],hfwd.src(2).vertno,'legacy'); 
          else
          [~,~,hrk] = intersect(hsrc(2).pinfo{sdiv.index},hfwd.src(2).vertno,'legacy');
          end
          Xls_plot(i).orient = hsrc(2).nn(hrk,:);
          %rr{i} = hsrc(2).rr(hrk,:);
      end
      nn(:,:,1) = Xls_plot(i).orient; nn = repmat(nn,[1,1,T]);
      Xls_plot(i).dip_or(:,:,1) = squeeze(nn(:,1,:)).*Xls_plot(i).dip;
      Xls_plot(i).dip_or(:,:,2) = squeeze(nn(:,2,:)).*Xls_plot(i).dip;
      Xls_plot(i).dip_or(:,:,3) = squeeze(nn(:,3,:)).*Xls_plot(i).dip;
      clear nn
    end
  else                                                                      % cortical case
    %% Backproject Cortical Scaled Estimates from Eigen Space to Physical Space
    patch_num = lk(c);                                                      % index of chosen patch
    cmode_inds = find(lk == patch_num);                                     % mode indices
    [~, ~, mode_num] = intersect(c, cmode_inds,'legacy');                   % index of chosen mode
    Vc = ico.patch(patch_num).V(:,mode_num);                                % right eigenvector for chosen mode - num_dipoles X 1
    Xls_cdip{i} = Vc*Xls_scaled(c,:);                                       % backproject modal estimates to dipoles
    Xls_plot(i).numdipoles = size(Vc,1);                                    % number of dipoles
    Xls_plot(i).cstrdip = currstr(c)/Xls_plot(i).numdipoles;                % current strength per dipole
    Xls_plot(i).mode = Xls_scaled(c,:);                                     % mode current estimate
    Xls_plot(i).dip = Xls_cdip{i};                                          % dipole current estimate
    Xls_plot(i).orient = cnn(ico.patch(patch_num).lk,:);                    % for orientation in MEG coordframe - use cfwd.source_nn(ico.patch(patch_num).lk,:) 
    nn(:,:,1) = Xls_plot(i).orient; nn = repmat(nn,[1,1,T]);
    Xls_plot(i).dip_or(:,:,1) = squeeze(nn(:,1,:)).*Xls_plot(i).dip;
    Xls_plot(i).dip_or(:,:,2) = squeeze(nn(:,2,:)).*Xls_plot(i).dip;
    Xls_plot(i).dip_or(:,:,3) = squeeze(nn(:,3,:)).*Xls_plot(i).dip;
    clear nn
  end
  
  %% Return Summary Time Courses for Plotting
  [Xls_plot(i).dip_res, Xls_plot(i).dip_resor] = vec_sum(Xls_plot(i).dip_or); % net dipole activity (amplitude and orientation) across all dipoles 
  Xls_plot(i).dip_resovern = Xls_plot(i).dip_res/Xls_plot(i).numdipoles;    % average dipole activity
  Xls_plot(i).dip_resnorm = Xls_plot(i).dip_res/(Xls_plot(i).cstrdip*Xls_plot(i).numdipoles); % average dipole activity (unitless)
end

%% Collate Time Series Across Modes for Same Patch or Subdivision
Xls_plot_modes = Xls_plot;
Xls_plot_modes = rmfield(Xls_plot_modes,{'mode','dip'});%,'dip_resor','dip_resovern','dip_resnorm'}); 
if sum(strimP) > 0
    strimP_reps = id_repmodes(reg_index(strimP(strimP>0),:), lk(strimP(strimP>0)), svolume, 's',strimP);
else
    strimP_reps = [];
end
if sum(ctrimP) > 0
    ctrimP_reps = id_repmodes(reg_index(ctrimP(ctrimP>0),:), lk(ctrimP(ctrimP>0)), cpatch, 'c',ctrimP);
else
    ctrimP_reps = [];
end
trimP_reps = [ctrimP_reps strimP_reps];

%% When multiple modes from same region, add up their contributions to dipole space
for i = 1:length(trimP_reps)
    trimind = trimP_reps{i}; 
    for j = 1:length(trimind)
        % note modes that are included
        ttrj = trimind(j);
        Xls_plot_modes(ttrj).mode_indices = trimind;
        Xls_plot_modes(ttrj).mode = [];

        % initialize dipole contributions
        Xls_plot_modes(ttrj).dip_or = zeros(size(Xls_plot_modes(ttrj).dip_or));
        Xls_plot_modes(ttrj).dip_res = zeros(size(Xls_plot_modes(ttrj).dip_res));
        
        % sum of contributions of all chosen modes to dipole space
        for kk = 1:length(trimind)
            Xls_plot_modes(ttrj).mode = [Xls_plot_modes(ttrj).mode; Xls_plot(trimind(kk)).mode];
            Xls_plot_modes(ttrj).dip_or = Xls_plot_modes(ttrj).dip_or + Xls_plot(trimind(kk)).dip_or;
        end
        
        % summary dipole activity (amplitude and orientation) across all dipoles
        [Xls_plot_modes(ttrj).dip_res, Xls_plot_modes(ttrj).dip_resor] = vec_sum(Xls_plot_modes(ttrj).dip_or);  %net
        Xls_plot_modes(ttrj).dip_resovern = Xls_plot_modes(ttrj).dip_res/Xls_plot_modes(ttrj).numdipoles;       % average dipole activity
        Xls_plot_modes(ttrj).dip_resnorm = Xls_plot_modes(ttrj).dip_res/(Xls_plot_modes(ttrj).cstrdip*Xls_plot_modes(ttrj).numdipoles); % average dipole activity (unitless)
    end
    clear trimind ttrj
end

end