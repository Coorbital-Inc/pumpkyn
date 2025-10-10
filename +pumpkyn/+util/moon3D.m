function [h,globe] = moon3D(posOffset,overlay,scale,h)
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

% Earth texture image
% Anything imread() will handle), but needs to be a 2:1 unprojected globe
    image_file = [fileparts(mfilename("fullpath")),filesep,'Moon.jpg'];

% Mean spherical earth
erad    = glb.EarthMoon.Rad; % equatorial radius (meters)
prad    = glb.EarthMoon.Rad; % polar radius (meters)
%erot    = 7.2921158553e-5; % earth rotation rate (radians/sec)

%% Create figure
if ~overlay
    figure('Color','k','Name','3D Earth Viewer','Renderer','opengl');
    hold on;
    % Turn off the normal axes
    set(gca, 'NextPlot','add','color','k');
    axis(gca,'off','equal','auto');
    % Set initial view
    view(0,30);
    axis equal;
    axis vis3d;
    h = gca;
end

if ~exist('h','var')
    h = gca;
end

%% Create wireframe globe
% Create a 3D meshgrid of the sphere points using the ellipsoid function
[x, y, z] = ellipsoid(0, 0, 0, erad, erad, prad, npanels);
x =  x.*scale + posOffset(1);
y =  y.*scale + posOffset(2);
z = -z.*scale + posOffset(3);
globe = surf(h, x, y, z, 'FaceColor', 'none', 'EdgeColor', 0.5*[1 1 1]); hold on;

%% Texturemap the globe
% Load Moon image for texture map

    cdata = imread(image_file);
     cmap = colormap(gray);
    cdata = pumpkyn.util.grs2rgb(cdata,cmap);
    % Set image as color data (cdata) property, and set face color to indicate
    % a texturemap, which Matlab expects to be in cdata. Turn off the mesh edges.
    set(globe, 'FaceColor', 'texturemap', 'CData', cdata, 'FaceAlpha', alpha, 'EdgeColor', 'none');
%    colormap gray;
globe.Annotation.LegendInformation.IconDisplayStyle = 'off';
end
