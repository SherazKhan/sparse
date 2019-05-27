function [leg_fig, clabel_fname, hlabel_fname] = ID_legends_labels(lk, trimP, clkmax, strimP, ctrimP, csrc, cfwd,...
         hsrc, hfwd, hipsurf_ind1, regname, ssubdiv, svolume, datapath, prefix, cdiv, sdiv_fwd_plot, hfwd_type, varargin)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Identify Anatomic Regions that Algorithm Converges On - Create Legends with these Regions
% Identify Label #s Corresponding to these Regions for Further Processing of Visuals
% Written by Pavitra Krishnaswamy (pavitrak@nmr.mgh.harvard.edu)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Initialize and Assign Inputs
if nargin == 19
    X = varargin{1};                                                        % only input X if simulated labels
end
leg_fig = cell(1, length(trimP));                                           % initialize MATLAB figure legends
init = zeros(1, length(trimP)); clblno = init; subcort_vol = init; clr = init;
cell_init = cell(1, length(trimP));  hlblno = cell_init; hlr = cell_init; 
 
%% Identify Region Names and Labels Corresponding to Estimates
for i = 1:length(trimP)
  s = strimP(i);                                                            % index wrt subcortical columns in joint srcspace
  c = ctrimP(i);                                                            % index wrt cortical columns in joint srscpace
  
  % Identify Cortical Legends and Labels
  if c > 0
      curr_patch = lk(c);                                                   % denote cortical patch
      if exist('X')                                                         % simulated labels 
          for j = 1:length(X)
           indc = curr_patch > min(X(j).lk) - 5 && curr_patch < max(X(j).lk) + 5;
           if indc == 1
              leg_fig{i}  = X(j).label;                                     % identify specific cortical region
           end
          end
      else                                                                  % no simulated labels
          leg_fig{i} = 'cortical';                                          % non-specific cortical region
      end
      cpatch_num = lk(c);                                                   % index of chosen patch
      [clblno(i), ~, ~, clr(i)] = patch_to_labloc(csrc,cfwd,cpatch_num);    % label # corresponding to patch
      %reg_index(i,:) = [cpatch_num, clr(i)];                               % unique identifier for every chosen patch
      
  % Identify Subcortical Legends and Labels
  elseif s > 0
      curr_patch = lk(s);                                                   % denote subcortical region
      for sreg = 1:length(regname)
          for sdiv = 1:length(ssubdiv{sreg})
              inds = find(svolume{sreg,sdiv} == curr_patch);                % ID anatomical region, subdivision index
              if isempty(inds) == 0
                  leg_fig{i} = strcat(regname{sreg}, '  subdiv', num2str(sdiv)); % display reg, subdiv
                  %reg_index(i,:) = [sreg, sdiv];                           % unique identifier for every chosen subdiv
                  slkmax = sum(lk<=clkmax);
                  if strfind(leg_fig{i},'lhipsurf')                         % left hippocampal patch
                      if isfield(sdiv_fwd_plot(s-slkmax), 'patchno')
                      hpatch_num = sdiv_fwd_plot(s-slkmax).patchno;         % index of chosen patches
                      else
                      hpatch_num = sdiv_fwd_plot(s-slkmax).index;           % index of chosen patch
                      end
                      for hhpp_len = 1:length(hpatch_num)
                      [hlblno{i}(hhpp_len), ~, ~, hlr{i}(hhpp_len)] = patch_to_labloc(hsrc,hfwd,hpatch_num(hhpp_len)); % label # for patch
                      end
                  elseif strfind(leg_fig{i},'rhipsurf')                     % right hippocampal patch
                      if isfield(sdiv_fwd_plot(s-slkmax), 'patchno')
                      hpatch_num = sdiv_fwd_plot(s-slkmax).patchno + length(hsrc(1).pinfo); % index of chosen patches
                      else
                      hpatch_num = sdiv_fwd_plot(s-slkmax).index + length(hsrc(1).pinfo);   % index of chosen patch
                      end
                      for hhpp_len = 1:length(hpatch_num)
                      [hlblno{i}(hhpp_len), ~, ~, hlr{i}(hhpp_len)] = patch_to_labloc(hsrc,hfwd,hpatch_num(hhpp_len)); % label # for patch
                      end
                  else
                      subcort_vol(i) = 1;                                   % subcortical volume
                  end
              end
          end
      end
  end  
  clear j indc sreg sdiv inds curr_patch 
end

%% Copy Relevant Labels to Current Folder
cnt = 1;
clabel_fname = cell(1,length(trimP)); hlabel_fname = clabel_fname;
clblno = squeeze(clblno); 
mkdir labels;
disp('Importing Relevant Label Files....');
for i = 1:length(trimP)
    if clblno(i) > 0
        clabel_fname{i} = strcat(num2str(clblno(i),'%06g'),'-',clr(i),'h.label');
        clabpath = strcat(datapath, 'label/',prefix,'-ico-',num2str(cdiv),'/',clr(i), 'h/');
        clabsource = strcat(clabpath,clabel_fname{i});        
        copyfile(clabsource,strcat(pwd,'/labels'));
        disp(strcat('Imported cortical label file .... ',clabel_fname{i}));
        cmgz_fname{i} = strcat(num2str(clblno(i),'%06g'),'-',clr(i),'h.mgz');
        unix_cmd{cnt} = ['mri_label2vol --label ', clabpath, clabel_fname{i}, ' --temp $SUBJECTS_DIR/$SUBJECT/mri/orig.mgz --identity --fillthresh 0 --o labels/', cmgz_fname{i},';'];
        cnt = cnt + 1;
    elseif ~isempty(hlblno{i})
        for hhpp_len = 1:length(hlblno{i})
        hlabel_fname{i,hhpp_len} = strcat(num2str(hlblno{i}(hhpp_len),'%06g'),'-',hlr{i}(hhpp_len),'h.mgz');
        hlabsource = strcat(datapath, 'label/',prefix,'_hip-',hfwd_type,'/',hlr{i}(hhpp_len), 'h/',hlabel_fname{i,hhpp_len});        
        copyfile(hlabsource,strcat(pwd,'/labels'));
        disp(strcat('Imported hippocampal label file .... ',hlabel_fname{i,hhpp_len}));
        end
    end
end

%% Display Labels and Message to User
disp('Subcortical volumes found...');
disp({leg_fig{subcort_vol==1}});
disp('Subcortical masks will be autogenerated....');
disp('Time to Generate Masks for Labels in Terminal....');
if sum(clblno) > 0
    disp(strcat('cortical labels......'));
    disp(clabel_fname(~cellfun('isempty',clabel_fname)));
else
    clabel_fname = [];
end
if sum(~cellfun(@isempty, hlblno)) > 0
    disp(strcat('hippocampal labels imported......'));
    disp(hlabel_fname(~cellfun('isempty',hlabel_fname)));
else
    hlabel_fname = [];
end
if sum(clblno) > 0
    disp(strcat('command to convert labels to masks ....'));
    disp(char(unix_cmd(~cellfun(@isempty,unix_cmd))));
    %disp('mri_label2vol --label #1.label --temp $SUBJECTS_DIR/$SUBJECT/mri/orig.mgz --identity --fillthresh 0 --o #1.mgz');
    %--fillthresh 0 gets every part of label, fillthresh 0.3 gets some part of it
    %--label-stat - to give values or colors to mask
end

end