function check_subdivs(prefix, datatype, meastype, ico_num, hfwd_type, append)

%% Paths and Prefixes
clc; close all;
dpath = '/autofs/eris/purdongp/users/pavitrak/sourceloc_data/';
datapath = strcat(dpath,datatype,'/',prefix);
addpath(genpath(datapath));
add_on_name = []; %add_on_name = '_bs-full'; 
load(strcat(datapath,'/fwd/',prefix,'_subc_ico',num2str(ico_num),'_',meastype,'_all_1g25m_SVDs',add_on_name,'.mat'),'sdiv_fwd'); 
whitened_svd = sdiv_fwd; clear sdiv_fwd;
load(strcat(datapath,'/fwd/',prefix,'_subc_ico',num2str(ico_num),'_',meastype,'_prewhite_volumes',add_on_name,'.mat'));
logFile  = [datapath, '/src/srcspaces_setto_ico',num2str(ico_num), '_details',add_on_name,'.txt'];          

%% Masks to Check Subdivisions in Regions of Interest - All Regions Processed
aseg = MRIread(strcat(datapath,'/mri/aseg.mgz'));
load(strcat(datapath,'/mri/T3.mat')); %T3 = [1 0 0 85.12; 0 1 0 128.00; 0 0 1 -128.00; 0 0 0 1.000];
v2ras1 = aseg.vox2ras1; mult = inv(v2ras1)*T3; 
sregion = {'lh','la','lp','lc','lt',...
           'rh','ra','rp','rc','rt',...
           'bs','lhipsurf','rhipsurf'};%'lic','ric',
if isempty(add_on_name)
    sregion{11} = 'bsred';
end
regcode = [17,18,12,11,10,16, 53,54,51,50,49,16,16, 17, 18];
newregcode = [1, 3, 6, 8, 12, 15, 17, 22, 28, 47, 37, 87, 14, 1, 3]; 
newregcode = repmat(newregcode, 1, 45);

for k = 1:length(sdiv_fwd)
    all_regions{k} = sdiv_fwd(k).reg_name;
end
c = 1; 
for i = 1:length(sregion)
    %% Volume ROI Masks
    %Initialize Region and Subdivisions
    regname{i} = sregion{i};                                                %name of ith region selected
    ind_sfwd = strncmp(regname{i}, all_regions, 4);
	ssubdiv{i} = find(ind_sfwd);                                        %subdivisions of region i
    ssubdiv{i}
     
    if isempty(strfind(regname{i},'surf'))
      % Initialize Mask
      m0{i} = aseg.vol==regcode(i); n = m0{i};                              %initialize aseg mask    
      
      for j = 1:length(ssubdiv{i})
        %Obtain MRI Coordinates - all are +/- sRAS + 129!
        rind = ssubdiv{i}(j)
        sRAS = sdiv_fwd(rind).srclocs(:,1:3)';                              %read sRAS coords region i subdiv j
        mricoords = mult(1:3,1:3)*sRAS;                                     %MRI voxels without additive factor
        yyn = mricoords(1,:) + mult(1,4);                                   %correcting for coordinate frame
        xxn = mricoords(2,:) + mult(2,4);                                   %correcting for coordinate frame
        zzn = mricoords(3,:) + mult(3,4);                                   %correcting for coordinate frame
        
        %Write Mask for Subdivision j Region i
        newn = n; newn(:,:,:) = 0;
        for mm = 1:length(xxn)
            newn(round(xxn(mm)), round(yyn(mm)), round(zzn(mm))) = newregcode(c);
        end
        c = c+1;
        A = aseg; A.vol = newn;
        MRIwrite(A,strcat(dpath,datatype,'/',prefix,'/mri/svolume_decomps/',regname{i},'_submask',num2str(j),'.mgz'));
        clear newn yyn xxn zzn sRAS vRAS mricoords A
      end
    else
     %% Surface ROI Masks
     %Initialize Generation of Masks from Surface ROIs in Label  Files
     hlr = sregion{i}(1:2);                                                 %left hippocampus or right hippocampus
     cnt = 1;
     dirName = strcat(datapath, '/label/',prefix,'_hip-',hfwd_type,'/',hlr, '/');%labels folder path
     destName = strcat(datapath, '/mri/svolume_decomps/');                  %MRI masks for subdivisions
     files = dir( fullfile(dirName,'*.label') );                            %# list all *.label files
     files = {files.name}';                                                 %# file names
    
     %General Purpose Masks
     unix_combined_label = 'mri_label2vol';
     for j = 1:numel(files)
      	hlabel_fname{j} = [dirName, files{j}];                              %# full path to source label file
        hmgz_fname{j} = [dirName,files{j}(1:end-6),'.mgz'];                 % full path to destination mask file          
        unix_cmd{cnt} = ['mri_label2vol --label ', hlabel_fname{j}, ' --temp $SUBJECTS_DIR/$SUBJECT/mri/orig.mgz --identity --fillthresh 0 --o ', hmgz_fname{j},';'];
        unix_combined_label = [unix_combined_label, ' --label ', hlabel_fname{j}];
        cnt = cnt + 1;
     end
     if ~exist(hmgz_fname{j},'file')
        disp(strcat('command to convert labels to masks ....'));
        disp([unix_combined_label ' --temp $SUBJECTS_DIR/$SUBJECT/mri/orig.mgz --identity --fillthresh 0 --o ', [datapath, '/mri/mask_', sregion{i},'.mgz'], ';']);
        disp('run the above in terminal'); pause;
        disp(char(unix_cmd(~cellfun(@isempty,unix_cmd))));
        disp('run the above in terminal'); pause;
        clear unix_cmd unix_combined_label;
     end
     
     %Masks Indicating Specific Subdivisions Made
     ind = strcmp(all_regions, sregion{i});
     indhip = find(ind, 1); num_patch = sum(ind);
     for ssdiv = 1:sum(ind)
        unix_cmd_fname{ssdiv} = 'mri_label2vol';
        if isfield(sdiv_fwd, 'patchno')
            patchind = sdiv_fwd(indhip+ssdiv-1).patchno;                    % patch grouping for hippocampus
        else
            patchind = sdiv_fwd(indhip+ssdiv-1).index;                      % no patch grouping for hippocampus
        end
        for pp = 1:length(patchind)
            unix_cmd_fname{ssdiv} = [unix_cmd_fname{ssdiv}, ' --label ', hlabel_fname{patchind(pp)}];
        end
     	hmgz_subdiv_fname{ssdiv} = [destName,sregion{i},'_submask',num2str(ssdiv),'.mgz'];% full path to destination mask file 
        unix_cmd_fname{ssdiv} = [unix_cmd_fname{ssdiv}, ' --temp $SUBJECTS_DIR/$SUBJECT/mri/orig.mgz --identity --fillthresh 0 --o ', hmgz_subdiv_fname{ssdiv},';'];     
     end  
     disp(char(unix_cmd_fname(~cellfun(@isempty,unix_cmd_fname))));
     disp('run the above in terminal'); pause;
     clear unix_cmd_fname 
    end
end

%% Summarize Anatomical and Electrophysiological Properties of Source Space - Only Regions Studied
efh = fopen(logFile,'w');                                                   % Log file to save computing times and diagnostics
fprintf(efh,'%s\t\t %s\t\t %s\t\t %s\t %s\t\t %s\t\t %s\t\t %s\n','Structure', 'SrcSpace Type', ...
        '#Dipoles','Size [mm3/2]','DMD [nAm/Size]', '#Sdivs', 'Size Per Div [mm3/2]','Modes per Div');
fclose(efh);
type_srcspace = repmat({'volume'}, 11, 1); type_srcspace{12} = 'surface'; type_srcspace{13} = 'surface';
load(strcat(datapath, '/fwd/curr_den.mat'),'subc_cd');
if strcmp(hfwd_type,'ico-1')
    hfwd_type = 'all';
end
for k = 1:length(whitened_svd)
    all_regions_white{k} = whitened_svd(k).reg_name;
end
for i = 1:length(sregion)
    ind = strcmp(all_regions_white, sregion{i}); 
    ssubdiv{i} = find(ind);                                            	    %subdivisions of region i
    if ~ isempty(ssubdiv{i})                            
      efh = fopen(logFile,'a'); 
      if strcmp(type_srcspace{i}, 'volume')                         
        load(strcat(datapath, '/src/', prefix, '_',sregion{i},'_geom.mat'));%geometry file
        reg_dims = round(volume);                                           %dimensions of region
        no_subdivs = length(ssubdiv{i});                                    %# subdivisions in region
        no_modes = median([whitened_svd(ind).numwhite_modes]);              %# modes per division in region (mean is not a round #)
        fprintf(efh,'%s\t\t\t\t %s\t\t\t\t %d\t\t\t %d\t\t\t\t %0.2f\t\t\t %d\t\t\t %0.2f\t\t\t\t %d\n',...
        sregion{i}, type_srcspace{i}, ndipoles, reg_dims, subc_cd(i)*(10^9), no_subdivs, reg_dims/no_subdivs, no_modes);
      else
        load(strcat(datapath, '/fwd/',prefix,'_',meastype,'_hipsurf_',hfwd_type,'_raw_SVDs',append,'.mat'));
        num_patch = length(ico.patch)/2;                                    %#patches in each hemisphere
        if strcmp(sregion{i}, 'lhipsurf')
            ndipoles = size(source.hemisph(1).rr, 1);   reg_dims = round(sum([ico.patch(1:num_patch).A])); %dimensions of region
        else
            ndipoles = size(source.hemisph(2).rr, 1);   reg_dims = round(sum([ico.patch(num_patch+1:end).A])); %dimensions of region
        end
        no_subdivs = length(ssubdiv{i});                                    %# subdivisions in region 
        no_modes = median([whitened_svd(ind).numwhite_modes]);              %# modes per division in region (mean is not a round #)
        fprintf(efh,'%s\t\t\t\t %s\t\t\t\t %d\t\t\t %d\t\t\t\t %0.2f\t\t\t %d\t\t\t %0.2f\t\t\t\t %d\n',...
        sregion{i}, type_srcspace{i}, ndipoles, reg_dims, subc_cd(i)*(10^9), no_subdivs, reg_dims/no_subdivs, no_modes);
      end
      fclose(efh);
    end
end

end