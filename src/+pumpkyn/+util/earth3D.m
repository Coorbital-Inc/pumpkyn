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
%                                                quality    'high' or
%                                                           'interactive'
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
        options.quality = 'high';
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
if ~isfield(options,'quality')
   options.quality = 'high';
end

quality = validatestring( ...
    options.quality, ...
    {'interactive','high'}, ...
    mfilename, ...
    'options.quality');
% if ~isfield(options,'logo')
%    options.logo = 'OD_LOGO.png'; 
% end
if ~isfield(options,'bgcolor') && ~strcmp(options.type,'BW')
   options.bgcolor = 'k'; 
elseif ~isfield(options,'bgcolor')
   options.bgcolor = 'k'; 
end
if ~strcmp(options.type,'BW') && strcmp(quality,'interactive')
    npanels = 180;
elseif ~strcmp(options.type,'BW')
    npanels = 500;   % Number of globe panels around the equator deg/panel = 360/npanels
else
    npanels = 360/20;
end
alpha   = 1; % globe transparency level, 1 = opaque, through 0 = invisible

% Earth texture image
if strcmp(options.type,'day')
    if strcmp(quality,'interactive') && ...
            isfile(fullfile( ...
            fileparts(mfilename('fullpath')), ...
            'earth-clouds-4k.jpg'))

        image_file = 'earth-clouds-4k.jpg';
    else
        image_file = 'earth-clouds-16k.jpg';
    end
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
    www = figure('Color',options.bgcolor,'NumberTitle','off','Name','3D Viewer'); 
end

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
[x,y,z] = getEarthMesh(erad,prad,npanels);
globe = surf(gca(www),options.posOffset(1)+x.*options.scale, ...
             options.posOffset(2)+y.*options.scale, ...
             options.posOffset(3)-z.*options.scale, ...
             'FaceColor', 'none', 'EdgeColor', 0.5*[1 1 1]); hold on;
shading(gca(www),'interp');

%% Create Atmosphere (Volumetric Onion Glow)
if options.atmos
    num_layers = 8;               
    max_alt_km = 100;             
    base_color = [0.3, 0.6, 1.0]; 
    
    rad_base = sqrt(x.^2 + y.^2 + z.^2);
    atmos = gobjects(num_layers, 1);
    
    for i = 1:num_layers
        fraction = (i / num_layers)^1.5; 
        layer_alt = max_alt_km * fraction;
        
        scale_ratio = (rad_base + layer_alt) ./ rad_base;
        
        % Apply user transformations exactly as done for the main globe
        x_atm = options.posOffset(1) + (x .* scale_ratio) .* options.scale;
        y_atm = options.posOffset(2) + (y .* scale_ratio) .* options.scale;
        z_atm = options.posOffset(3) - (z .* scale_ratio) .* options.scale;
        
        layer_alpha = 0.05 * (1 - fraction)^2; 
        
        atmos(i) = surf(gca(www), x_atm, y_atm, z_atm, ...
            'FaceColor', base_color, ...
            'FaceAlpha', layer_alpha, ...
            'EdgeColor', 'none'); hold on;
    end
end

%% Texturemap the globe
% Load Earth image for texture map
if ~strcmp(options.type,'none')
    cdata = readEarthTexture( ...
        image_file,strcmp(quality,'interactive'));
    if strcmp(options.type,'infrared')
       idx = sum(cdata,3) >= 240*3;
       for td=1:size(cdata,3)
               tmp = cdata(:,:,td);
          tmp(idx) = 0;
          cdata(:,:,td) = tmp;
       end
    end
    
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
    haPos = get(gca,'position');
    ha2=axes('position',[1-haPos(1),haPos(2),0.1,0.1,]);
    image(imread(options.logo));
    axis(ha2,'equal');
    set(ha2,'handlevisibility','off','visible','off');
end

%% Insert Star Background:
if options.stars
    bh = star3D(www);
    uistack(bh, 'bottom');
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
    else
        hLight = lightangle(gca(www),180-20,0);
    end
    
    material(gca(www),'shiny');
    
    % --- Solid Body (Globe) ---
    globe.FaceLighting = 'gouraud';     
    globe.SpecularStrength = 0.05;      
    globe.SpecularExponent = 10;        
    globe.DiffuseStrength = 1.0;        
    globe.AmbientStrength = 0.05;       
    globe.BackFaceLighting = 'unlit';   
    globe.EdgeLighting = 'none';
    globe.Annotation.LegendInformation.IconDisplayStyle = 'off';
    
    % --- Atmosphere (Volumetric Shells) ---
    if options.atmos
        for i = 1:length(atmos)
            atmos(i).FaceLighting = 'gouraud';     
            atmos(i).SpecularStrength = 0.0;       
            atmos(i).DiffuseStrength = 1.0;        
            atmos(i).AmbientStrength = 0.4;        
            atmos(i).BackFaceLighting = 'lit';     
            atmos(i).EdgeLighting = 'none';
            atmos(i).Annotation.LegendInformation.IconDisplayStyle = 'off';
        end
    end
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

function [x,y,z] = getEarthMesh(erad,prad,npanels)
%% Cache the most recently used base Earth mesh.

persistent cachedEquatorialRadius cachedPolarRadius cachedPanelCount
persistent cachedX cachedY cachedZ

cacheMatches = ...
    ~isempty(cachedPanelCount) && ...
    cachedPanelCount == npanels && ...
    cachedEquatorialRadius == erad && ...
    cachedPolarRadius == prad;

if ~cacheMatches
    [cachedX,cachedY,cachedZ] = ellipsoid( ...
        0,0,0,erad,erad,prad,npanels);

    cachedEquatorialRadius = erad;
    cachedPolarRadius = prad;
    cachedPanelCount = npanels;
end

x = cachedX;
y = cachedY;
z = cachedZ;

end

function cdata = readEarthTexture(imageFile,useCache)
%% Cache only interactive textures to avoid retaining a 16K image in RAM.

persistent cachedFile cachedImage

if useCache && ...
        ~isempty(cachedFile) && ...
        strcmp(cachedFile,imageFile)

    cdata = cachedImage;
    return;
end

cdata = imread(imageFile);

if useCache
    cachedFile = imageFile;
    cachedImage = cdata;
end

end
