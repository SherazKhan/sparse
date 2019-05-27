%% Paths and File Names
clear; clc; close all;
addpath /usr/local/freesurfer/stable5_3_0/matlab;
dpath = '/autofs/eris/purdongp/users/pavitrak/sourceloc_data/AEP/';
subject = 'cat_004';
cd(strcat(dpath,subject,'/mri'));
mask_names = {'lic','ric'};
col_masks = [22 22]; ds =1;
m0 = cell(size(mask_names));
saveon = 1;

for i = 1:length(m0)
    %% Read Mask Corresponding to a Label
    fname = strcat(mask_names{i},'.mgz');
    reg = MRIread(fname);
    m0{i} = reg.vol==1;
    n = m0{i}; 

    %% Read Voxel Indices Corresponding to Mask
    f = find(n == 1);
    [x,y,z]=ind2sub(size(n),f); %gives indices where n has 1's

    %% Approximate Mask = Subsample of Actual Mask
    newn = n; newn(:,:,:)= 0;
    xxn = x(1:ds:end); yyn = y(1:ds:end); zzn = z(1:ds:end);
    for mm = 1:length(xxn)
        newn(xxn(mm), yyn(mm), zzn(mm)) = 1; %col_masks(i)
    end
	    
    %% Visualize Masks in MATLAB
    figure(i), set(gcf,'color','white'); 
    plot3(y,x,z,'k.'); hold on; 
    plot3(yyn, xxn, zzn,'r*'); hold on;
    xlabel('Pos-Ant Voxels'); ylabel('Left-Right Voxels'); zlabel('Inf-Sup Voxels');
    title(mask_names{i}); axis equal; axis vis3d; 
    if mod(i,2)
        view([-145 10]);
    end

    %% Write Mask to Visualize in Freeview
    A = reg; B = reg;
    A.vol = n; B.vol = newn; 
    if saveon
        %MRIwrite(A,strcat('mask_',mask_names{i},'_orig.mgz'));
        MRIwrite(B,strcat('mask_',mask_names{i},'.mgz'));
    end

    %% For Source Space -- Transform from Voxels to Surface RAS Coordinates
    vRAS = reg.vox2ras1*[yyn'; xxn'; zzn'; ones(size(yyn'))];
    %vox2ras1 accounts for matlab voxel numbering which starts from 1
    %terminal - "mne_collect_transforms --mgh T1.mgz" in terminal and obtain T1
    load('T3.mat'); %T3 = [1 0 0 -0.5; 0 1 0 -0.00; 0 0 1 -0.00; 0 0 0 1.000];
    sRAS = T3\vRAS; %inv(T3)*vRAS

    %% Write Source Space File - All Three Orientations
    coord = sRAS(1:3,:)'; cnt = 1;
    for kk = 1:size(coord,1)
        coord_or(cnt,:) = [coord(kk,:) 1 0 0];
        coord_or(cnt+1,:) = [coord(kk,:) 0 1 0];
        coord_or(cnt+2,:) = [coord(kk,:) 0 0 1];
        cnt = cnt + 3; 
    end
    savename = strcat(dpath,subject,'/src/',subject,'_',mask_names{i},'.txt');
    if saveon
        save(savename, 'coord_or', '-ASCII'); %type(savename);
    end

    %% Save Geometries and Subdivisions
    ndipoles = size(coord,1);           %not accounting orientation
    volume = ndipoles;
    %ndipoles_or = size(coord_or,1);     %includes orientations
    if saveon
        save(strcat(dpath, subject, '/src/',subject,'_',mask_names{i},'_geom.mat'),'volume','ndipoles');
    end
    clear n x y z xxn yyn zzn newn B A vRAS sRAS coord coord_or cnt volume ndipoles %ndipoles_or;
end