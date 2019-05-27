function [sregname, sdiv_fwd, best_raw_sdiv_modes, best_white_sdiv_modes] ...
    = svolume_decomp_sdiv(subjid, meastype, datapath, volumes_fname, megCw, eegCw, megsel, eegsel, nthresh, p_sing)
%% Load Files
load(volumes_fname,'sregname', 'sdiv_fwd');
load(strcat(datapath, '/fwd/curr_den.mat'),'subc_cd');                      %current densities from literature

%% Pick out Good Channels and Whiten Fwd Solutions - by Subdivision
for k = 1:length(sdiv_fwd)
    if strcmp(meastype,'meg')                                               %MEG only case
        sdiv_fwd(k).whiteG = megCw*sdiv_fwd(k).rawG(megsel,:);              %whitened forward solution
        sdiv_fwd(k).gdchnums = megsel;                                      %good channel numbers
        %sdiv_fwd(k).gdchnames = megch;
    elseif strcmp(meastype,'eeg')                                           %EEG only case
        sdiv_fwd(k).whiteG = eegCw*sdiv_fwd(k).rawG(eegsel,:);              %whitened forward solution
        sdiv_fwd(k).gdchnums = eegsel;                                      %good channel numbers
        %sdiv_fwd(k).gdchnames = eegch;
    else                                                                    %MEG-EEG combined
        sdiv_fwd(k).whiteG = cat(1, megCw*sdiv_fwd(k).rawG(megsel,:), eegCw*sdiv_fwd(k).rawG(eegsel,:));
        sdiv_fwd(k).gdchnums = [megsel, eegsel];                            %good channel numbers
        %sdiv_fwd(k).gdchnames = [megch eegch];
    end
end

%% If IC Exists: Merge the Left and Right IC Before SVD
for k = 1:length(sdiv_fwd)
    ind(k) = isempty(strfind(sdiv_fwd(k).reg_name,'ic'));                   %find ric and lic
end
ic_ind = find(ind==0);                                                      %the indices for IC
if ~isempty(ic_ind)
    sdiv_noic = sdiv_fwd(ind==1);                                           %sdiv with no IC
    sdiv_fwd = [sdiv_noic sdiv_fwd(ind==0)];                                %sdiv with IC wherein both IC terms are together
    sdiv_fwd(end-1).reg_name = 'ic'; sdiv_fwd(end-1).index = 1;             %rename and reindex
    sdiv_fwd(end).reg_name = 'ic';   sdiv_fwd(end).index = 2;               %rename and reindex
    ric_name = strcmp(sregname, 'ric'); sregname = sregname(~ric_name);     %no need to consider ric separate from lic
    sregname{strcmp(sregname,'lic')} = 'ic';                                %rename to IC
    ic_name = strcmp(sregname,'ic'); sregname = [sregname(~ic_name) sregname(ic_name)];
    subc_cd = subc_cd(~ric_name); subc_cd = [subc_cd(~ic_name) subc_cd(ic_name)];%current densities reordered
end

%% Analysis by Subdivision
for k = 1:length(sdiv_fwd)
  disp(strcat('Computing Volume Decomposition for...   ', sdiv_fwd(k).reg_name, '...  index... ', num2str(k)));
  
  % Reduced Order SVDs of Fwd Solutions 
  ploton = 0;  savemem = 0;                                                 %parameters for reduced SVD
  [Gr,Ur,Sr,Vr,pr,nmr_r{k},~] = nmra(sdiv_fwd(k).rawG,ploton,savemem,p_sing,nthresh);   %svds of raw forward solution
  [Gw,Uw,Sw,Vw,pw,nmr_w{k},~] = nmra(sdiv_fwd(k).whiteG,ploton,savemem,p_sing,nthresh); %svds of raw forward solution
  sdiv_fwd(k).raw_nmra = [nmr_r{k}];                                        %mean rep. accuracy of raw fwd solution
  sdiv_fwd(k).numraw_modes = pr;                                            %raw case - #modes needed for 95% NMRA
  sdiv_fwd(k).white_nmra = [nmr_w{k}];                                      %mean rep. accuracy of white fwd solution
  sdiv_fwd(k).numwhite_modes = pw;                                          %white case - #modes needed for 95% NMRA
  
  % Current Strength
  regnum = strcmp(sdiv_fwd(k).reg_name, sregname);                          %index of fwd region 
  sdiv_fwd(k).curr_den = subc_cd(regnum);                                   %current density of this region
  if isempty(sdiv_fwd(k).area)                                              %scale current density by volume
      curr_str(k) = subc_cd(regnum)*sdiv_fwd(k).volume;                     %current strength for all modes used (no power loss)
      %curr_str(k) = subc_cd(regnum)*sdiv_fwd(k).volume/sdiv_fwd(k).ndipole;%current strength per dipole
  elseif isempty(sdiv_fwd(k).volume)                                        %scale current density by surface area
      curr_str(k) = subc_cd(regnum)*sdiv_fwd(k).area;                       %current strength for all modes used (no power loss)
      %curr_str(k) = subc_cd(regnum)*sdiv_fwd(k).area/sdiv_fwd(k).ndipole;  %current strength per dipole
  end
  sdiv_fwd(k).currstr = curr_str(k);                                        %strengths of subdiv for this region
  
  % Store Appropriately Scaled Eigenmodes
  sdiv_fwd(k).rawGstheta = curr_str(k)*Gr;                                  %raw case - implicitly scaled by S, scale by currstr
  sdiv_fwd(k).rawUs = Ur;                                                   %raw case - left singular vector - sensor modes
  sdiv_fwd(k).rawS = Sr;                                                    %raw case - singular values
  sdiv_fwd(k).rawVs = Vr;                                                   %raw case - right singular vector - spatial modes
  sdiv_fwd(k).whiteGstheta = curr_str(k)*Gw;                                %white case - implicitly scaled by S, scale by currstr
  sdiv_fwd(k).whiteUs = Uw;                                                 %white case - left singular vector - sensor modes
  sdiv_fwd(k).whiteS = Sw;                                                  %white case - singular values
  sdiv_fwd(k).whiteVs = Vw;                                                 %white case - right singular vector - spatial modes
end

%% Plot NMRAs Across Subdivisions and Choose Effective #Modes/Subdivision
for ii = 1:p_sing 
    % NMRA Across Region and Median Energy Captured by ii Modes
    for k = 1:length(sdiv_fwd)
        if length(nmr_r{k}) >= ii
        nmra_r(k) = nmr_r{k}(ii);                                           %iith NMRA for each region
        end
        if length(nmr_w{k}) >= ii
        nmra_w(k) = nmr_w{k}(ii);                                           %iith NMRA for each region
        end
    end
    
    % Plot NMRA for Raw G
    if exist('nmra_r')
    Er(ii) = median(nmra_r);                                                %Median of iith NMRA across regions
    figure(3), set(gcf, 'color', 'white') 
    subplot(2,round(p_sing/2),ii); hist(nmra_r)
    grid on; title({'subcortical subdivisions',strcat(' p = ',num2str(ii),' modes')})
    xlabel('NMRA'); ylabel('Number of Volumes')
    hold on; plot(Er(ii),0,'*r')
    end
    
    % Plot NMRA for White G
    if exist('nmra_w')
    Ew(ii) = median(nmra_w);                                                %Median of iith NMRA across regions
    figure(4), set(gcf, 'color', 'white') 
    subplot(2,round(p_sing/2),ii); hist(nmra_w)
    grid on; title({'subcortical subdivisions',strcat(' p = ',num2str(ii),' modes')})
    xlabel('NMRA'); ylabel('Number of Volumes')
    hold on; plot(Ew(ii),0,'*r')
    end
    
    clear nmra_r nmra_w
end
figure(3), h = gcf; print(h,'-dpdf',strcat(subjid, '_raw_nmras_sdiv_',meastype,'_scort'));
figure(4), h = gcf; print(h,'-dpdf',strcat(subjid, '_white_nmras_sdiv_',meastype,'_scort'));
best_raw_sdiv_modes = find(round(Er*100) >= nthresh,1);                     %for same #modes across regions
best_white_sdiv_modes = find(round(Ew*100) >= nthresh,1);                   %for same #modes across regions
close all;

% sdiv_fwd had regname, subdiv index, rawG, src locs/inds, volumes, and ndipoles fields
% sdiv_fwd now has whiteG, good channel numbers/names, raw/white Gstheta/Us/S/p/Vs, curr_den and curr_str fields 
end