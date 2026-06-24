function [h,globe,hLight] = moon3D(posOffset,overlay,scale,h)
%% Purpose:
% Provides an interface to view the 3D moon in space.  This allows for
% custom visualization of satellite orbits around the earth. You can
% overlay a track onto the 3D plot by plotting the ECEF position of a
% satellite in km.
%
%% Inputs:
%
%% Revision History:
%  Darin Koblick        Modified to support additional earth marble images
%  Copyright 2025 Coorbital, Inc.
%% -------------------- Begin Code Sequence -------------------------------
glb = pumpkyn.util.getConst();

if nargin == 0
   posOffset = [0,0,0];
   overlay = false;
   scale = 1;
end

   npanels = 100;   % Number of globe panels around the equator deg/panel = 360/npanels
   alpha   = 1.0;   % globe transparency level, 1 = opaque, through 0 = invisible
bump_scale = 1;     % Exaggeration factor for the topography. 1 = True to life, 10 = Visible

% Earth texture image
    image_file = [fileparts(mfilename("fullpath")),filesep,'Moon_4k.jpg'];
    disp_file  = [fileparts(mfilename("fullpath")),filesep,'Moon_4k_disp.tif'];

% Mean spherical earth
erad    = glb.EarthMoon.Rad; % equatorial radius (meters)
prad    = glb.EarthMoon.Rad; % polar radius (meters)

%erot    = 7.2921158553e-5; % earth rotation rate (radians/sec)

%% Create figure
if ~overlay
    figure('Color','k','Name','3D Moon Viewer');
    hold on;
    % Turn off the normal axes
    set(gca, 'NextPlot','add','color','k');
    axis(gca,'off');
    % Set initial view
    view(0,30);
    h = gca;
end

if ~exist('h','var')
    h = gca;
end

% Force 1:1:1 aspect ratio and freeze it for 3D rotation (CRITICAL FIX)
axis(h, 'equal');
axis(h, 'vis3d');

%% Create wireframe globe
% Create a 3D meshgrid of the sphere points using the ellipsoid function
[x, y, z] = ellipsoid(0, 0, 0, erad, erad, prad, npanels);

%% Apply Surface Displacement ---
if exist(disp_file, 'file')
    % 1. Read the floating-point TIFF (values are in km)
    disp_map = double(imread(disp_file)); 
    % 2. Resize the displacement map to perfectly match the ellipsoid vertex grid
    disp_resized = imresize(disp_map, size(x));
    % 3. Convert displacement from km to meters, and apply our visual exaggeration
    disp_meters = disp_resized .* 1000 .* bump_scale;
    % 4. Calculate unit vectors for every point on the sphere (pointing outward)
    rad_orig = sqrt(x.^2 + y.^2 + z.^2);
    x_norm = x ./ rad_orig;
    y_norm = y ./ rad_orig;
    z_norm = z ./ rad_orig;
    % 5. Push the vertices out along their normal vectors by the displacement amount
    rad_new = rad_orig + disp_meters./1000;
    x = x_norm .* rad_new;
    y = y_norm .* rad_new;
    z = z_norm .* rad_new;
else
    warning('Displacement map TIFF not found. Rendering a perfectly smooth sphere.');
end


%% Apply positional scaling and offsets
    x = -x.*scale + posOffset(1);
    y = -y.*scale + posOffset(2);
    z = -z.*scale + posOffset(3);
globe = surf(h, x, y, z, 'FaceColor', 'none', 'EdgeColor', 0.5*[1 1 1]); hold on;


%% Texturemap the globe
% Load Moon image for texture map

    cdata = imread(image_file);
    if size(cdata, 3) == 1
         cmap = colormap(gray);
        cdata = pumpkyn.util.grs2rgb(cdata, cmap);
    end

    % Set image as color data (cdata) property, and set face color to indicate
     set(globe, 'FaceColor', 'texturemap', 'CData', cdata, 'FaceAlpha', alpha, 'EdgeColor', 'none');
     globe.Annotation.LegendInformation.IconDisplayStyle = 'off';
     set(gca,'clipping','off');

% Optional: Change shading to make the 3D bumps interact with figure lighting
 lighting gouraud;
 hLight = light('Position', [1 0 0], 'Style', 'infinite', 'Color', [1 1 1]);
 material([0.5, 0.8, 0.0]);

end
