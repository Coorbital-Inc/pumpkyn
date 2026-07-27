function [h,globe,hLight] = moon3D( ...
    posOffset,overlay,scale,h,varargin)
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

parser = inputParser;
parser.FunctionName = 'pumpkyn.util.moon3D';

addParameter( ...
    parser,'Quality','high', ...
    @(value) (ischar(value) && isrow(value)) || ...
    (isstring(value) && isscalar(value)));

addParameter( ...
    parser,'UseDisplacement',[], ...
    @(value) isempty(value) || ...
    (islogical(value) && isscalar(value)) || ...
    (isnumeric(value) && isscalar(value) && ismember(value,[0 1])));

addParameter( ...
    parser,'AddLighting',true, ...
    @(value) (islogical(value) && isscalar(value)) || ...
    (isnumeric(value) && isscalar(value) && ismember(value,[0 1])));

parse(parser,varargin{:});

quality = validatestring( ...
    parser.Results.Quality, ...
    {'interactive','high'}, ...
    parser.FunctionName, ...
    'Quality');

useDisplacement = parser.Results.UseDisplacement;

if isempty(useDisplacement)
    useDisplacement = strcmp(quality,'high');
else
    useDisplacement = logical(useDisplacement);
end

addLighting = logical(parser.Results.AddLighting);

if strcmp(quality,'interactive')
    npanels = 72;
else
    npanels = 100;
end

   alpha   = 1.0;   % globe transparency level, 1 = opaque, through 0 = invisible
bump_scale = 1;     % Exaggeration factor for the topography. 1 = True to life, 10 = Visible

% Earth texture image
    assetPath = fileparts(mfilename('fullpath'));

    if strcmp(quality,'interactive') && ...
            isfile(fullfile(assetPath,'Moon_2k.jpg'))

        image_file = fullfile(assetPath,'Moon_2k.jpg');
    else
        image_file = fullfile(assetPath,'Moon_4k.jpg');
    end

    disp_file = fullfile(assetPath,'Moon_4k_disp.tif');

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
[x,y,z] = getMoonMesh( ...
    erad,prad,npanels, ...
    useDisplacement,disp_file,bump_scale);


%% Apply positional scaling and offsets
    x = -x.*scale + posOffset(1);
    y = -y.*scale + posOffset(2);
    z = -z.*scale + posOffset(3);
globe = surf( ...
    h,x,y,z, ...
    'FaceColor','none', ...
    'EdgeColor',0.5*[1 1 1]);

hold(h,'on');


%% Texturemap the globe
% Load Moon image for texture map

    cdata = readMoonTexture( ...
        image_file,strcmp(quality,'interactive'));
    if size(cdata, 3) == 1
         cmap = colormap(gray);
        cdata = pumpkyn.util.grs2rgb(cdata, cmap);
    end

    % Set image as color data (cdata) property, and set face color to indicate
     set(globe, 'FaceColor', 'texturemap', 'CData', cdata, 'FaceAlpha', alpha, 'EdgeColor', 'none');
     globe.Annotation.LegendInformation.IconDisplayStyle = 'off';
     set(h,'clipping','off');

% Configure only the Moon surface instead of changing material properties
% for every surface already present in the axes.
hLight = gobjects(0);

if addLighting
    set(globe, ...
        'FaceLighting','gouraud', ...
        'AmbientStrength',0.5, ...
        'DiffuseStrength',0.8, ...
        'SpecularStrength',0.0, ...
        'BackFaceLighting','unlit');

    hLight = light( ...
        h, ...
        'Position',[1 0 0], ...
        'Style','infinite', ...
        'Color',[1 1 1]);
else
    globe.FaceLighting = 'none';
end

end

function [x,y,z] = getMoonMesh( ...
    erad,prad,npanels,useDisplacement,dispFile,bumpScale)
%% Cache the most recently requested lunar surface geometry.

persistent cachedKey cachedX cachedY cachedZ

cacheKey = sprintf( ...
    '%.15g|%.15g|%d|%d|%.15g', ...
    erad,prad,npanels,useDisplacement,bumpScale);

if ~isempty(cachedKey) && strcmp(cachedKey,cacheKey)
    x = cachedX;
    y = cachedY;
    z = cachedZ;
    return;
end

[x,y,z] = ellipsoid( ...
    0,0,0,erad,erad,prad,npanels);

if useDisplacement
    if isfile(dispFile)
        displacementMap = double(imread(dispFile));
        displacement = imresize(displacementMap,size(x));

        originalRadius = sqrt(x.^2+y.^2+z.^2);
        newRadius = originalRadius+displacement*bumpScale;

        x = x./originalRadius.*newRadius;
        y = y./originalRadius.*newRadius;
        z = z./originalRadius.*newRadius;
    else
        warning( ...
            'pumpkyn:moon3D:MissingDisplacementMap', ...
            ['Displacement map TIFF not found. Rendering a ', ...
             'perfectly smooth sphere.']);
    end
end

cachedKey = cacheKey;
cachedX = x;
cachedY = y;
cachedZ = z;

end

function cdata = readMoonTexture(imageFile,useCache)
%% Cache only the reduced interactive texture.

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
