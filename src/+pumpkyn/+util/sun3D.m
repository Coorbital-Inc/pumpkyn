function [h,globe] = sun3D(posOffset,overlay,scale,h)
%% Purpose:
% Provides an interface to view the 3D Sun in space.  This allows for
% custom visualization of satellite orbits around the earth. You can
% overlay a track onto the 3D plot by plotting the ECEF position of a
% satellite in km.
%
%% Inputs:
%
%% Revision History:
%  Darin Koblick        
%  Copyright 2026 Coorbital, Inc.
%% -------------------- Begin Code Sequence -------------------------------
glb = pumpkyn.util.getConst();

if nargin == 0
   posOffset = [0,0,0];
   overlay = false;
   scale = 1;
end

npanels = 32;    % A solid Sun needs only a modest mesh for animation
sunColor = [1.0,0.72,0.05];
% Mean spherical earth
erad    = glb.Sun.Rad; % equatorial radius (meters)
prad    = glb.Sun.Rad; % polar radius (meters)
%erot    = 7.2921158553e-5; % earth rotation rate (radians/sec)

%% Create figure
if ~overlay
    figure('Color','k','Name','3D Sun Viewer');
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

%% Create Sun globe
% Create a 3D meshgrid of the sphere points using the ellipsoid function
[x, y, z] = ellipsoid(0, 0, 0, erad, erad, prad, npanels);
x =  x.*scale + posOffset(1);
y =  y.*scale + posOffset(2);
z = -z.*scale + posOffset(3);
globe = surf( ...
    h,x,y,z, ...
    'FaceColor',sunColor, ...
    'FaceLighting','none', ...
    'FaceAlpha',1, ...
    'EdgeColor','none');
globe.Annotation.LegendInformation.IconDisplayStyle = 'off';
end
