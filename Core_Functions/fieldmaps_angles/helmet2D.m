function helmet2D(proj_approach, helmet_fname, evoked_fname)

%% Create Helmet File
mnepath = '/autofs/eris/purdongp/users/pavitrak/matlab/mne/share/mne/mne_analyze/';
s = mne_read_bem_surfaces(strcat(mnepath,'306m.fif'),true);
Vertices3D = s.rr;                          % 3D vertices where sensors located
Faces3D = s.tris;                           % Faces based on Even Triang in 3D
save(helmet_fname,'Vertices3D','Faces3D');  % Save vertices and faces

%% Project 3D to 2D
switch proj_approach
    % Cartesian to Sphere Transform
    case 'cart_sph'
        x = Vertices3D(:,1);                % 3D cartesian coordinate 1
        y = Vertices3D(:,2);                % 3D cartesian coordinate 2
        z = Vertices3D(:,3);                % 3D cartesian coordinate 3
        z = z - max(z);                     % Take out Z Direction
        [TH,PHI,R] = cart2sph(x, y, z);     % 3D cartesian to spherical coordinates
        PHI(PHI < 0.001) = 0.001;           % Remove the too smal values for PHI
        R2 = R./cos(PHI).^.2;               % Flat projection - stretch out for 3rd dim
        [X,Y] = pol2cart(TH,R2);            % 2D polar to cartesian coordinates
        Vertices2D =[X Y];                  % 2D vertices

    % Multidimensional Scaling
    case 'mdimsc'
        D = pdist(s.rr,'Mahalanobis');      % Distances between all pairs of sensors
        Vertices2D = mdscale(D,2);          % Scale to maintain distances between sensors
    otherwise 
        disp('Unknown Method');             % Error
end
Faces2D = delaunay(Vertices2D);             % Faces based on Even Triang in 2D
save(helmet_fname,'Vertices2D','Faces2D','-append');

%% Write Layout File
%% Write Layout File
nChan = size(Vertices2D,1);                 % Number of channels
slno = 1:nChan;                             % Serial No.
info = fiff_read_meas_info(evoked_fname);   % Measurement Info        
chnames = info.ch_names(1:nChan);           % Channel Names
M = cell(nChan, 6);
for i = 1:nChan
    M{i,1} = num2str(slno(i),'%03d');
    M{i,2} = Vertices2D(i,1)*200;
    M{i,3} = Vertices2D(i,2)*170;
    M{i,4} = 4.00;
    M{i,5} = 3.00;
    M{i,6} = chnames{i};
end
fileID = fopen(strcat(mnepath, 'magnesWH3600.lout'),'w');
formatSpec = '%s\t %.2f\t %.2f\t %.2f\t %.2f\t %s\n';
[nrows,ncols] = size(M);
for row = 1:nrows
    fprintf(fileID,formatSpec,M{row,:});
end
fclose(fileID);
type(strcat(mnepath, 'lout/magnesWH3600.lout'));

end