function auto_s_srcspace(subject, datatype)

%% General Paths and Parameters
clc; close all;
addpath /usr/local/freesurfer/stable5_3_0/matlab;
dpath = strcat('/autofs/eris/purdongp/users/pavitrak/sourceloc_data/',datatype,'/');
cpath = '/autofs/eris/purdongp/users/pavitrak/sourceloc_scripts/';
datapath = strcat(dpath,subject);
addpath(strcat(cpath,subject,'-code/srcfwd_final'));
saveon = 1; ploton = 1;

%% Region Code from Freeview Label System
fullname = {...
    'Left-Thalamus-Proper','Right-Thalamus-Proper','Left-Caudate','Right-Caudate',...
    'Left-Putamen','Right-Putamen','Left-Hippocampus','Right-Hippocampus', ...
    'Left-Amygdala','Right-Amygdala','Brain-Stem'};
tregion = {...
    'lthalamus_prop','rthalamus_prop','lcaudate','rcaudate',...
    'lputamen','rputamen', 'lhippo','rhippo','lamygdala', 'ramygdala','brainstem'};
load(strcat(cpath,'ScSPIGH/source_spaces/sregcode.mat'), 'region', 'reg_code','code');
for i = 1:length(tregion)
    s(i) = strmatch(tregion{i},region);
    tcode{i} = code{s(i)};
end

%% Read ASEG and Initialize
aseg = MRIread(strcat(datapath, '/mri/aseg.mgz')); 
load(strcat(datapath,'/mri/T3.mat')); 	%T3 = [1 0 0 85.12; 0 1 0 128.00; 0 0 1 -128.00; 0 0 0 1.000];
v2ras1 = aseg.vox2ras1; mult = v2ras1\T3;                                   %RAS to voxel conversion 
mask_names = {'lt','rt','lc','rc','lp','rp','lh','rh','la','ra','bs'};
col_masks = [10 10 11 11 12 12 17 17 18 18 16]; ds = 1;
m0 = cell(size(mask_names));
rng('default'); num_bs = 4;                                                 %for reduction of brainstem srcspace

for i = 1:length(tregion)
    %% Read Mask Corresponding to a Label
    m0{i} = aseg.vol==tcode{i};
    n = m0{i}; 
    
    %% Read Voxel Indices Corresponding to Mask
    f = find(n == 1);
    [x,y,z] = ind2sub(size(n),f);                                           %gives indices where n has 1's
    
    %% Approximate Mask = Subsample of Actual Mask
    newn = n; newn(:,:,:)= 0; 
    xxn = x(1:ds:end); yyn = y(1:ds:end); zzn = z(1:ds:end);
    for mm = 1:length(xxn)
        newn(xxn(mm), yyn(mm), zzn(mm)) = 1;                                %col_masks(i);
    end
    
    %% Visualize Masks in MATLAB
    figure(i), set(gcf,'color','white'); 
    plot3(y,x,z,'k.'); hold on; 
    plot3(yyn, xxn, zzn,'r*'); hold on;
    xlabel('Pos-Ant Voxels'); ylabel('Left-Right Voxels'); zlabel('Inf-Sup Voxels');
    title(tregion{i}); axis equal; axis vis3d; 
    if mod(i,2)
        view([-145 10]);
    end
   
    %% Write Mask to Visualize in Freeview
    A = aseg; B = aseg;
    A.vol = n; B.vol = newn;
    if saveon
        %MRIwrite(A,strcat(datapath,'/mri/mask_',mask_names{i},'_orig.mgz'));
        MRIwrite(B,strcat(datapath,'/mri/mask_',mask_names{i},'.mgz'));
        %save(strcat(datapath,'/mri/','check_contin',num2str(i),'.mat'),'yyn','xxn','zzn','newn');
    end
    
    %% For Source Space -- Transform from Voxels to Surface RAS Coordinates
    vRAS = aseg.vox2ras1*[yyn'; xxn'; zzn'; ones(size(yyn'))];
    %vox2ras1 accounts for matlab voxel numbering which starts from 1
    %terminal - "mne_collect_transforms --mgh T1.mgz" in terminal and obtain T1
    sRAS = T3\vRAS; %inv(T3)*vRAS
    
    %% Write Source Space File - All Three Orientations
    coord = sRAS(1:3,:)'; cnt = 1;
    for kk = 1:size(coord,1)
        coord_or(cnt,:) = [coord(kk,:) 1 0 0];
        coord_or(cnt+1,:) = [coord(kk,:) 0 1 0];
        coord_or(cnt+2,:) = [coord(kk,:) 0 0 1];
        cnt = cnt + 3; 
    end
    savename = strcat(datapath,'/src/',subject,'_',mask_names{i},'.txt');
    if saveon
        save(savename, 'coord_or', '-ASCII');                               %type(savename);
    end
    
    %% Save Geometries and Subdivisions
    %terminal - "asegstats2table --i $SUBJECTS_DIR/$SUBJECT/stats/aseg.stats --transpose --tablefile aseg_stats.txt"
    mask_vols = importdata(strcat(dpath,subject,'/stats/aseg_stats.txt'));
    volnum = find(strncmp(fullname{i}, mask_vols.rowheaders, 15)==1);
    volume = mask_vols.data(volnum);
    ndipoles = size(coord,1);                                               %not accounting orientation
    if saveon
        save(strcat(datapath, '/src/',subject,'_',mask_names{i},'_geom.mat'),'volume','ndipoles');
    end
    
    %% Brainstem Reduction: Src Space and Mask
    if strcmp(tregion{i}, 'brainstem')
      % Reduce Src Space - sRAS Coordinates
      srclocs_bsfull = importdata(strcat(datapath,'/src/',subject,'_',mask_names{i},'.txt'));       %Source Locations
      [srclocs, ndipoles, volume] = reduce_brainstem(datapath, subject, srclocs_bsfull, ndipoles, volume);
      %[srclocs, ndipoles, volume, pick_inds] = reduce_brainstem(srclocs_bsfull, ndipoles, ...
      %volume, num_bs, strcat(datapath, '/src/',mask_names{i},'-reduction'), ploton);
      save(strcat(datapath, '/src/',subject,'_',mask_names{i},'red.txt'), 'srclocs', '-ASCII');     %Reduced source locations
      save(strcat(datapath, '/src/',subject,'_',mask_names{i},'red_geom.mat'),'ndipoles','volume'); %Reduced source volume and # dipoles

      % Reduced SrcSpace - Obtain MRI Coordinates 
      sRAS_red = srclocs(1:3:end,1:3)';                                     %MRI coordinates (sRAS)
      mricoords_red = mult(1:3,1:3)*sRAS_red;                               %MRI voxels without additive factor
      yyn_red = mricoords_red(1,:) + mult(1,4);                             %correcting for coordinate frame
      xxn_red = mricoords_red(2,:) + mult(2,4);                             %correcting for coordinate frame
      zzn_red = mricoords_red(3,:) + mult(3,4);                             %correcting for coordinate frame
    
      % Reduced SrcSpace - MRI Mask 
      newn_red = n; newn_red(:,:,:) = 0; %newn_red = zeros(size(aseg.vol)); %initialize mask
      for mm = 1:length(xxn_red)
        newn_red(round(xxn_red(mm)), round(yyn_red(mm)), round(zzn_red(mm))) = 1;
      end
      Ared = aseg; Ared.vol = newn_red;
      MRIwrite(Ared,strcat(datapath,'/mri/mask_',mask_names{i},'red.mgz'));
      clear newn_red yyn_red xxn_red zzn_red sRAS_red mricoords_red Ared
    end
    clear n x y z xxn yyn zzn newn vRAS  sRAS coord coord_or cnt mask_vols volume ndipoles;
end

%% Run SrcFwd Commands
disp(['cd /autofs/eris/purdongp/users/pavitrak/sourceloc_scripts/',subject,'-code/srcfwd_final/SrcFwd_Shell_Scripts']);
disp('script ~/scort_srcfwd_log.txt'); 
disp('./scort_srcfwd_exec');
disp('exit');
disp('script ~/cort_srcfwd_log.txt'); 
disp('./cort_srcfwd_exec');
disp('exit');
disp('script ~/hip_srcfwd_log.txt'); 
disp('./hip_srcfwd_exec');
disp('exit');
disp(['after it runs: type Exit and Copy the exec and log files to', ...
     subject,'-code/fwd_final/SrcFwd_Shell_Scripts and sourceloc_data/',datatype,'/',subject,'/fwd/term_logs']);
 pause;

end