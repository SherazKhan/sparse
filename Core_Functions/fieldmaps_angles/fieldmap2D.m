function cax = fieldmap2D(helmet_fname, proj_approach, evoked_fname, data, nContours, cax)

%% Create/Load 2D Helmet File
if exist(helmet_fname,'file') ~= 2
    helmet2D(proj_approach, helmet_fname, evoked_fname)
end
load(helmet_fname,'Faces3D', 'Vertices2D');%'Vertices3D', 'Faces2D'

%% Initialize Figure
figure, set(gcf,'color','white');
hAxes = gca;                            % Axes for Contour/NoseEars Plotting
e_alpha = 0;%0.2;                       % Edge Thickness
f_alpha = 1;                            % Face Transparency
e_color = 'none';%'k';                  % Edge Color
s_marker = 'ro';                        % Sensor Location Marker

%% Basic Patches
hold on; 
patch_obj = patch('Faces', Faces3D, 'Vertices', Vertices2D, 'EdgeColor', e_color, 'LineStyle', 'none', 'EdgeAlpha', e_alpha, ...
'FaceVertexCData', data, 'FaceColor', 'interp', 'BackfaceLighting', 'lit','facealpha',f_alpha);

%% Contours (Based on Brainstorm)
if nContours > 0
    hold on;
    tricontour(Vertices2D, Faces3D, data, nContours, hAxes); 
end

%% Nose Ears Background
hold on;
radii = [Vertices2D(:,2); Vertices2D(:,1)];
PlotNoseEars(hAxes, (max(radii)-min(radii))/4, 1);

%% Sensor Locations
%hold on;
%scatter(Vertices2D(:,1),Vertices2D(:,2), s_marker);

%% Sensor Names
%for i=1:size(Vertices2D,1)
%    text(Vertices2D(i,1), Vertices2D(i,2),0,num2str(i));
%end

%% Axes Settings, Scales and View Orientations
colormap('jet'); colorbar; 
set(gca,'fontsize', 12, 'fontweight', 'b');
if isempty(cax)
    cax = [-max(abs(data)) max(abs(data))];
end
set(gca,'clim',[cax(1) cax(2)]);
axis('equal'); 
axis('off');
view(0,89.5);

end