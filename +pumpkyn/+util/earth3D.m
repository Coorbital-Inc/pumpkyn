function [www,globe,hLight,cdata] = earth3D(options,www)
%% Purpose:
% Provides an interface to view the 3D earth in space.  This allows for
% custom visualization of satellite orbits around the earth. You can
% overlay a track onto the 3D plot by plotting the ECEF position of a
% satellite in km.
%
%% Inputs:
%  options              struct                      type:
%
%                                                   'clouds'    Earth image
%                                                               with cloud
%                                                               cover
%
%                                                   'night'     Earth image
%                                                               with city
%                                                               lights
%
%                                                   'day'       Earth image
%                                                               during
%                                                               sunlight
%
%                                           'some clouds'       Earth image
%                                                               with part
%                                                               cloud cover
%                                                              
%                                                   animate:
%                                                  true         animate the
%                                                               rotation of
%                                                               the earth
%                                                               
%                                                   false       do not
%                                                               animate
%                                                               earth's
%                                                               rotation
%
%                                                   bgcolor     'k'
%                                                               'w'
%
%
%                                                      logo:
%                                                           if specified,
%                                                           this will
%                                                           display a PNG
%                                                           logo.
%
%                                               posOffset   If field is
%                                                           specified then
%                                                           the Earth will
%                                                           be offset by
%                                                           this amount in
%                                                           ECI space
%                                                           [x,y,z]
%
%                                                scale      Scaling ratio
%                                                           of
%                                                           Earth Radius
%                                               
%
%% Revision History:
%  Darin Koblick        Modified to support additional earth marble images
%  Darin Koblick        Modifed Additional earth marble images for speed
%  Darin Koblick        Added julian date field for lighting effects
%  Copyright 2025 Coorbital, Inc.
%% -------------------- Begin Code Sequence -------------------------------

if ~exist('options','var')
         options.type = 'day'; %'BW'; %'day';
      options.animate = false; %true;
      options.bgcolor = 'k'; %'w'; %'k';
        options.stars = false; %false; %true;
    options.posOffset = [0,0,0];
        options.scale = 1;
        options.AddShading = true;
        options.overlay = false;
end
if ~isfield(options,'animate')
     options.animate = false; 
end
if ~isfield(options,'type')
        options.type = 'night'; 
end
if ~isfield(options,'stars')
       options.stars = false; 
end
if ~isfield(options,'posOffset')
   options.posOffset = [0 0 0];
end
if ~isfield(options,'scale')
       options.scale = 1; 
end
if ~isfield(options,'AddShading')
   options.AddShading = true; 
end
if ~isfield(options,'overlay')
   options.overlay = false; 
end
if ~isfield(options,'WGS84')
   options.WGS84 = false; 
end
if ~isfield(options,'atmos')
   options.atmos = true; 
end
% if ~isfield(options,'logo')
%    options.logo = 'OD_LOGO.png'; 
% end
if ~isfield(options,'bgcolor') && ~strcmp(options.type,'BW')
   options.bgcolor = 'k'; 
elseif ~isfield(options,'bgcolor')
   options.bgcolor = 'k'; 
end

if ~strcmp(options.type,'BW')
    npanels = 500;   % Number of globe panels around the equator deg/panel = 360/npanels
else
    npanels = 360/20;
end
alpha   = 1; % globe transparency level, 1 = opaque, through 0 = invisible

% Earth texture image
% Anything imread() will handle), but needs to be a 2:1 unprojected globe
if strcmp(options.type,'day')
    image_file = 'world.topo.200412.3x5400x2700.png'; %'BMNG_world.topo.bathy.200405.3.2048x1024.jpg';
elseif strcmp(options.type,'night')
     image_file = 'land_lights_16384.tif';
elseif strcmp(options.type,'clouds')
    image_file = 'cloud_combined_2048.tif';  
elseif strcmp(options.type,'some clouds')
    image_file = 'land_ocean_ice_cloud_2048.tif';
elseif strcmp(options.type,'BW')
    image_file = 'outline-black-white-world-political.jpg';
elseif strcmp(options.type,'infrared')
    image_file = 'InfraredEarth.png';
elseif strcmp(options.type,'altEarth')
    image_file = 'alternateEarth.bmp';
end

%% Put full path into image file:
image_file = [fileparts(mfilename('fullpath')),filesep,image_file];

% Mean spherical earth
if options.WGS84
    erad = 6378.1370;
    prad = 6356.75231424500;
else
     glb = pumpkyn.util.getConst();
    erad = glb.Earth.Rad; % equatorial radius (km)
    prad = glb.Earth.Rad; % polar radius (km)
end

%% Create figure if not supplied
if nargin <2
    www = figure('Color',options.bgcolor,'NumberTitle','off','Name','3D Viewer', ...
                 'Renderer','opengl'); 
end

%set(www,'Renderer','opengl')
%opengl hardware;
% Turn off the normal axes
if ~options.overlay
hold on;
set(gca(www), 'NextPlot','add', 'Visible','off');
% Set initial view
view(gca(www),0,30);
axis(gca(www),'equal');
set(gca(www),'clipping','off');
end
%% Create wireframe globe
% Create a 3D meshgrid of the sphere points using the ellipsoid function
  [x,y,z] = ellipsoid(0, 0, 0, erad, erad, prad, npanels);
    globe = surf(gca(www),options.posOffset(1)+x.*options.scale, ...
                 options.posOffset(2)+y.*options.scale, ...
                 options.posOffset(3)-z.*options.scale, ...
                 'FaceColor', 'none', 'EdgeColor', 0.5*[1 1 1]); hold on;
    shading(gca(www),'interp');

%% Create Atmosphere
if options.atmos
    [x, y, z] = ellipsoid(0, 0, 0, erad+100, erad+100, prad+100, npanels);
    atmos = surf(gca(www),options.posOffset(1)+x.*options.scale, ...
                 options.posOffset(2)+y.*options.scale, ...
                 options.posOffset(3)-z.*options.scale, ...
                 'FaceColor', [0.7 0.7 1], 'EdgeColor', 'none','facealpha',0.2); hold on;
end
%% Texturemap the globe
% Load Earth image for texture map
if ~strcmp(options.type,'none')
    cdata = imread(image_file);
    if strcmp(options.type,'infrared')
       idx = sum(cdata,3) >= 240*3;
       for td=1:size(cdata,3)
               tmp = cdata(:,:,td);
          tmp(idx) = 0;
          cdata(:,:,td) = tmp;
       end

    end
    % Set image as color data (cdata) property, and set face color to indicate
    % a texturemap, which Matlab expects to be in cdata. Turn off the mesh edges.
    if ~strcmp(options.type,'BW')
        set(globe, 'FaceColor', 'texturemap', 'CData', cdata, 'FaceAlpha', alpha, 'EdgeColor', 'none');
    else
        set(globe, 'FaceColor', 'texturemap', 'CData', cdata, 'FaceAlpha', alpha, ...
            'EdgeColor', [0.2 0.2 0.2],'LineStyle','--');
    end
end

%% Insert a logo over the 3D plot
if isfield(options,'logo')
    uistack(gca,'bottom');
    % Creating a new axes for the logo on the current axes
    % To create the logo at the bottom left corner of the plot use
    % the next two lines
    haPos = get(gca,'position');
    ha2=axes('position',[1-haPos(1),haPos(2),0.1,0.1,]);
    % To place the logo at the bottom left corner of the figure window
    % uncomment the line below and comment the above two lines
    % ha2=axes('position',[0, 0, .1,.04,]);
    % Adding a LOGO to the new axes
    % The logo file(jpeg, png, etc.) must be placed in the working path
    image(imread(options.logo));
    % Turn the handlevisibility off so that we don't inadvertently plot
    % into the axes again. Also, make the axes invisible
    axis(ha2,'equal');
    set(ha2,'handlevisibility','off','visible','off');
end

%% Insert Star Background:

if options.stars
    %ah = axes('unit', 'normalized', 'position', [0 0 1 1]);
    % import the background image and show it on the axes
    %bg = imread('starmap_g8k.jpg'); 
    %sbg = imagesc(bg);
    %set(sbg, 'AlphaData', 0.75);
    % prevent plotting over the background and turn the axis off
    %set(ah,'handlevisibility','off','visible','off')
    % making sure the background is behind all the other uicontrols
    %uistack(ah, 'bottom');
    %Map the 3D star scene:
    bh = star3D(www);
    uistack(bh, 'bottom');
    %Set the camera position:
    set(gca,'cameraviewanglemode','manual');
    axis vis3d;
    axis equal;
    zoom(100);
end
hLight = [];
if ~strcmp(options.type,'BW') && options.AddShading
    %Adjust 3D lighting:
    if isfield(options,'jd0')
        %Determine the sun position:
    [rSun,vSun] = sunPosVel(options.jd0);
        %Convert to ECEF:
           rSun = ECItoECEF(options.jd0,rSun,vSun,vSun.*0,2);
        %Find the apparent azimuth and elevation angle:
        [az,el] = cart2sph(rSun(1),rSun(2),rSun(3));
         hLight = lightangle(gca(www),270-az*180/pi,el*180/pi);
        %Put in an arrow pointing toward the sun:
        %quiver3(0,0,0,8000*rSun(:,1)./vmag(rSun,2), ...
        %              8000*rSun(:,2)./vmag(rSun,2), ...
        %              8000*rSun(:,3)./vmag(rSun,2));
    else
        hLight = lightangle(gca(www),180-20,20);
    end
    material(gca(www),'shiny');
    globe.FaceLighting = 'flat';
    globe.SpecularStrength = 0.1;
    globe.SpecularExponent = 0.2;
    globe.SpecularColorReflectance = 0.4;
    globe.DiffuseStrength = 1;
    globe.AmbientStrength = 0.25;
    globe.BackFaceLighting = 'unlit';
    globe.EdgeLighting = 'flat';
    atmos.FaceLighting = 'flat';
    atmos.SpecularStrength = 0.5;
    atmos.SpecularExponent = 0.1;
    atmos.SpecularColorReflectance = 0.2;
    atmos.DiffuseStrength = 1;
    atmos.AmbientStrength = 0.1;
    atmos.BackFaceLighting = 'lit';
    atmos.EdgeLighting = 'flat';
    atmos.Annotation.LegendInformation.IconDisplayStyle = 'off';
    globe.Annotation.LegendInformation.IconDisplayStyle = 'off';
end

if options.animate
    while true
            for i=0:359
                if ishandle(globe)
                    pause(0.01);
                    view(i,30);
                    drawnow;
                else
                    return;
                end
            end
    end
end

end
