function hStarSphere = stars3D( ...
    fig,ax,backgroundMode,starImageFile,jd0)
%% Purpose:
%
% Display an equirectangular starfield on the inside of a celestial sphere.
%
% The starfield is drawn in a separate background axes. Foreground camera
% changes are copied to the starfield using one-way listeners. Changes to
% the starfield cannot modify the foreground camera.
%
%% Inputs:
%
% fig               handle    Target figure
%
% ax                handle    Foreground three-dimensional axes
%
% backgroundMode    char      'stars' or 'black'
%
% starImageFile     char      Equirectangular 2:1 starfield image
%
% jd0               double    Julian date used to orient the starfield
%
%% Output:
%
% hStarSphere       handle    Starfield surface handle
%

%% Self test

if nargin == 0
          jd0 = pumpkyn.util.juliandate( ...
              '03/21/2000 00:00:00');

        lStar = 389703.264829278;
        muStar = 0.012150585609624;

    [fig,globe] = pumpkyn.cr3bp.showEarth( ...
        jd0,lStar,muStar);

              ax = ancestor(globe,'axes');

    starImageFile = fullfile( ...
        fileparts(mfilename('fullpath')), ...
        'starmap_16k.jpg');

    if ~isfile(starImageFile)
        starImageFile = '';
    end

    view(ax,45,20);
    axis(ax,'equal');
    axis(ax,'vis3d');

    camtarget(ax,[-muStar,0,0]);
    camtarget(ax,'manual');

    hStarSphere = pumpkyn.util.stars3D( ...
        fig,ax,'stars',starImageFile,jd0);

    rotate3d(fig,'on');
    set(fig,'CurrentAxes',ax);

    drawnow;
    return;
end

%% Defaults

if nargin < 3 || isempty(backgroundMode)
    backgroundMode = 'stars';
end

if nargin < 4
    starImageFile = '';
end

if nargin < 5 || isempty(jd0)
    jd0 = 2451545.0;
end

%% User-adjustable starfield settings

starViewAngle = 65;    % Starfield field of view in degrees
sphereScale = 100;      % Sphere radius relative to visible scene size
spherePadding = 1.01;  % Padding applied to star-axes limits

numLongitude = 721;    % Sphere mesh longitude samples
numLatitude = 361;     % Sphere mesh latitude samples

fallbackWidth = 2048;
fallbackHeight = 1024;

%% Validate inputs

if ~isgraphics(fig,'figure')
    error('fig must be a valid MATLAB figure handle.');
end

if ~isgraphics(ax,'axes')
    error('ax must be a valid MATLAB axes handle.');
end

if ancestor(ax,'figure') ~= fig
    error('ax must belong to the supplied figure.');
end

if ~strcmpi(backgroundMode,'stars') && ...
        ~strcmpi(backgroundMode,'black')

    error( ...
        'backgroundMode must be either ''stars'' or ''black''.');
end

%% Preserve the foreground before changing any graphics objects

foreground = captureForegroundAxes(ax);

% The cleanup restores the original foreground even if stars3D encounters
% an error while constructing the starfield.
foregroundCleanup = onCleanup(@() ...
    restoreForegroundAxes(fig,ax,foreground));

%% Remove any existing star background

% This occurs after capturing the foreground because an older background
% may contain bidirectional links created by a previous implementation.
removeExistingBackground(fig);

%% Select the background mode

if strcmpi(backgroundMode,'black')
    % Restore the foreground before selecting the opaque black background.
    clear foregroundCleanup;

    set(fig,'Color','k');
    set(ax,'Color','k');

    set(fig,'CurrentAxes',ax);

    hStarSphere = gobjects(0);
    return;
end

set(fig,'Color','k');

%% Load the starfield image

starImageFile = resolveImagePath(starImageFile);

if ~isempty(starImageFile) && isfile(starImageFile)
    starImage = readStarImage(starImageFile);
else
    warning( ...
        ['The requested starfield image was not found. ', ...
         'Using a generated starfield instead.']);

    starImage = generateStarfield( ...
        fallbackWidth,fallbackHeight);
end

validateImageProjection(starImage);

%% Determine the celestial-sphere size

cameraDistance = norm( ...
    foreground.CameraPosition - ...
    foreground.CameraTarget);

sceneSpan = norm([ ...
    diff(foreground.XLim), ...
    diff(foreground.YLim), ...
    diff(foreground.ZLim)]);

referenceDistance = max([ ...
    cameraDistance, ...
    sceneSpan, ...
    eps]);

sphereCenter = foreground.CameraTarget;
sphereRadius = sphereScale*referenceDistance;

%% Construct the celestial sphere

[xSphere,ySphere,zSphere] = createSphereGeometry( ...
    jd0, ...
    sphereCenter, ...
    sphereRadius, ...
    numLongitude, ...
    numLatitude);

%% Create the background axes

hStarAxes = axes( ...
    'Parent',fig, ...
    'Units',foreground.Units, ...
    'Position',foreground.Position, ...
    'Color','k', ...
    'Visible','off', ...
    'HitTest','off', ...
    'PickableParts','none', ...
    'HandleVisibility','off', ...
    'Projection','perspective', ...
    'DataAspectRatio',[1 1 1], ...
    'DataAspectRatioMode','manual', ...
    'SortMethod','depth', ...
    'Clipping','off', ...
    'Tag','stars3DBackgroundAxes');

uistack(hStarAxes,'bottom');
uistack(ax,'top');

% The foreground axes must be transparent so the background axes shows.
set(ax,'Color','none');

%% Configure the starfield camera

axesRadius = spherePadding*sphereRadius;

set(hStarAxes, ...
    'XLim',sphereCenter(1) + axesRadius*[-1 1], ...
    'YLim',sphereCenter(2) + axesRadius*[-1 1], ...
    'ZLim',sphereCenter(3) + axesRadius*[-1 1], ...
    'XLimMode','manual', ...
    'YLimMode','manual', ...
    'ZLimMode','manual', ...
    'CameraPosition',foreground.CameraPosition, ...
    'CameraTarget',foreground.CameraTarget, ...
    'CameraUpVector',foreground.CameraUpVector, ...
    'CameraViewAngle',starViewAngle, ...
    'CameraPositionMode','manual', ...
    'CameraTargetMode','manual', ...
    'CameraUpVectorMode','manual', ...
    'CameraViewAngleMode','manual', ...
    'Projection','perspective', ...
    'DataAspectRatio',[1 1 1], ...
    'DataAspectRatioMode','manual', ...
    'SortMethod','depth', ...
    'Clipping','off');

axis(hStarAxes,'vis3d');
hold(hStarAxes,'on');

%% Draw the starfield

hStarSphere = surf( ...
    hStarAxes, ...
    xSphere, ...
    ySphere, ...
    zSphere, ...
    starImage, ...
    'FaceColor','texturemap', ...
    'EdgeColor','none', ...
    'FaceAlpha',1, ...
    'FaceLighting','none', ...
    'BackFaceLighting','unlit', ...
    'Clipping','off', ...
    'HitTest','off', ...
    'PickableParts','none', ...
    'HandleVisibility','off', ...
    'Tag','stars3DBackground');

hStarSphere.Annotation.LegendInformation.IconDisplayStyle = 'off';

%% Add one-way foreground-to-starfield synchronization

% These listeners observe only the foreground axes. Starfield changes
% therefore cannot propagate back and alter view(ax,az,el).
synchronizeStarAxes(ax,hStarAxes);

hCameraListeners = { ...
    addlistener(ax,'CameraPosition','PostSet', ...
        @(~,~) synchronizeStarAxes(ax,hStarAxes)); ...
    addlistener(ax,'CameraTarget','PostSet', ...
        @(~,~) synchronizeStarAxes(ax,hStarAxes)); ...
    addlistener(ax,'CameraUpVector','PostSet', ...
        @(~,~) synchronizeStarAxes(ax,hStarAxes)); ...
    addlistener(ax,'Position','PostSet', ...
        @(~,~) synchronizeStarAxes(ax,hStarAxes))};

% Listener objects must remain in memory for synchronization to continue.
hStarSphere.UserData = struct( ...
    'StarAxes',hStarAxes, ...
    'CameraListeners',{hCameraListeners}, ...
    'StarViewAngle',starViewAngle);

%% Restore the exact foreground view

% Clearing the cleanup object invokes restoreForegroundAxes.
clear foregroundCleanup;

drawnow;

end

function foreground = captureForegroundAxes(ax)
%% Purpose:
%
% Capture the foreground view before creating the starfield.
%

foreground.CameraPosition = ax.CameraPosition;
foreground.CameraTarget = ax.CameraTarget;
foreground.CameraUpVector = ax.CameraUpVector;
foreground.CameraViewAngle = ax.CameraViewAngle;
foreground.Projection = ax.Projection;

foreground.XLim = ax.XLim;
foreground.YLim = ax.YLim;
foreground.ZLim = ax.ZLim;

foreground.DataAspectRatio = ...
    ax.DataAspectRatio;

foreground.Units = ax.Units;
foreground.Position = ax.Position;

end

function restoreForegroundAxes(fig,ax,foreground)
%% Purpose:
%
% Restore and freeze the foreground view exactly as it was supplied.
%

if ~isgraphics(ax,'axes')
    return;
end

set(ax, ...
    'CameraPosition',foreground.CameraPosition, ...
    'CameraTarget',foreground.CameraTarget, ...
    'CameraUpVector',foreground.CameraUpVector, ...
    'CameraViewAngle',foreground.CameraViewAngle, ...
    'Projection',foreground.Projection, ...
    'XLim',foreground.XLim, ...
    'YLim',foreground.YLim, ...
    'ZLim',foreground.ZLim, ...
    'DataAspectRatio',foreground.DataAspectRatio, ...
    'CameraPositionMode','manual', ...
    'CameraTargetMode','manual', ...
    'CameraUpVectorMode','manual', ...
    'CameraViewAngleMode','manual', ...
    'XLimMode','manual', ...
    'YLimMode','manual', ...
    'ZLimMode','manual', ...
    'DataAspectRatioMode','manual', ...
    'Color','none');

if isgraphics(fig,'figure')
    set(fig,'CurrentAxes',ax);
end

end

function synchronizeStarAxes(ax,hStarAxes)
%% Purpose:
%
% Copy foreground orientation and position to the starfield axes.
%

if ~isgraphics(ax,'axes') || ...
        ~isgraphics(hStarAxes,'axes')
    return;
end

set(hStarAxes, ...
    'CameraPosition',ax.CameraPosition, ...
    'CameraTarget',ax.CameraTarget, ...
    'CameraUpVector',ax.CameraUpVector, ...
    'Position',ax.Position);

end

function removeExistingBackground(fig)
%% Purpose:
%
% Delete any previously created stars3D background.
%

oldStarAxes = findall( ...
    fig, ...
    'Type','axes', ...
    'Tag','stars3DBackgroundAxes');

if ~isempty(oldStarAxes)
    validAxes = oldStarAxes(isgraphics(oldStarAxes));

    if ~isempty(validAxes)
        delete(validAxes);
    end
end

% Remove surfaces created by earlier single-axes implementations.
oldStarSphere = findall( ...
    fig, ...
    'Type','surface', ...
    'Tag','stars3DBackground');

if ~isempty(oldStarSphere)
    validSurfaces = ...
        oldStarSphere(isgraphics(oldStarSphere));

    if ~isempty(validSurfaces)
        delete(validSurfaces);
    end
end

end

function starImageFile = resolveImagePath(starImageFile)
%% Purpose:
%
% Resolve a relative image path against the directory containing stars3D.
%

if isempty(starImageFile) || isfile(starImageFile)
    return;
end

mFilePath = fileparts(mfilename('fullpath'));

localImageFile = fullfile( ...
    mFilePath,starImageFile);

if isfile(localImageFile)
    starImageFile = localImageFile;
end

end

function starImage = readStarImage(starImageFile)
%% Purpose:
%
% Read the starfield image and convert it to RGB uint8 data.
%

[starImage,colorMap] = imread(starImageFile);

if ~isempty(colorMap)
    starImage = uint8( ...
        255*ind2rgb(starImage,colorMap));
end

if ndims(starImage) == 2
    starImage = repmat(starImage,[1 1 3]);

elseif size(starImage,3) > 3
    starImage = starImage(:,:,1:3);
end

if isa(starImage,'uint16')
    starImage = uint8( ...
        double(starImage)./65535.*255);

elseif isa(starImage,'double') || ...
        isa(starImage,'single')

    if max(starImage(:)) <= 1
        starImage = uint8(255*starImage);
    else
        starImage = uint8(starImage);
    end
end

end

function validateImageProjection(starImage)
%% Purpose:
%
% Check that the image is approximately a 2:1 equirectangular map.
%

imageAspectRatio = ...
    size(starImage,2)/size(starImage,1);

if abs(imageAspectRatio-2) > 0.15
    warning( ...
        ['The starfield image does not have an approximately 2:1 ', ...
         'aspect ratio. It may not use an equirectangular ', ...
         'projection and could appear distorted.']);
end

end

function [xSphere,ySphere,zSphere] = createSphereGeometry( ...
    jd0,sphereCenter,sphereRadius,numLongitude,numLatitude)
%% Purpose:
%
% Construct the celestial sphere in the Earth-Moon CR3BP frame.
%

%% Construct unit directions in Galactic coordinates

longitude = linspace( ...
    pi,-pi,numLongitude);

latitude = linspace( ...
    pi/2,-pi/2,numLatitude);

[longitudeGrid,latitudeGrid] = ...
    meshgrid(longitude,latitude);

cosLatitude = cos(latitudeGrid);

rGalactic = [ ...
    (cosLatitude(:).*cos(longitudeGrid(:))).'; ...
    (cosLatitude(:).*sin(longitudeGrid(:))).'; ...
     sin(latitudeGrid(:)).'];

%% Rotate Galactic coordinates into J2000 coordinates

CG2I = [ ...
   -0.0548755604,  0.4941094279, -0.8676661490; ...
   -0.8734370902, -0.4448296300, -0.1980763734; ...
   -0.4838350155,  0.7469822445,  0.4559837762];

rJ2000 = CG2I*rGalactic;

%% Rotate J2000 coordinates into the CR3BP frame

[r12,v12] = pumpkyn.util.planetPosVel( ...
    jd0,'Earth','Moon');

r12 = reshape(r12(1:3),1,3);
v12 = reshape(v12(1:3),1,3);

xHat = r12./norm(r12);

zHat = cross(r12,v12);
zHat = zHat./norm(zHat);

yHat = cross(zHat,xHat);
yHat = yHat./norm(yHat);

CR2I = [ ...
    xHat(:), ...
    yHat(:), ...
    zHat(:)];

rCR3BP = CR2I.'*rJ2000;

%% Scale and position the sphere

sphereSize = size(longitudeGrid);

xSphere = sphereCenter(1) + ...
    sphereRadius*reshape( ...
    rCR3BP(1,:),sphereSize);

ySphere = sphereCenter(2) + ...
    sphereRadius*reshape( ...
    rCR3BP(2,:),sphereSize);

zSphere = sphereCenter(3) + ...
    sphereRadius*reshape( ...
    rCR3BP(3,:),sphereSize);

end

function starImage = generateStarfield(imageWidth,imageHeight)
%% Purpose:
%
% Generate a deterministic fallback equirectangular starfield.
%

starImage = zeros( ...
    imageHeight,imageWidth,3,'uint8');

randomStream = RandStream( ...
    'mt19937ar','Seed',271828);

numStars = 6000;

starX = randi( ...
    randomStream,imageWidth,numStars,1);

starY = randi( ...
    randomStream,imageHeight,numStars,1);

brightness = ...
    rand(randomStream,numStars,1).^4;

brightness = uint8( ...
    30 + 225*brightness);

for kk = 1:numStars
    xx = starX(kk);
    yy = starY(kk);

    starColor = ...
        double(brightness(kk))*[1.00 0.97 0.92];

    starImage(yy,xx,:) = uint8(reshape( ...
        min(255,starColor),1,1,3));

    if brightness(kk) > 180
        haloColor = uint8( ...
            0.25*starColor);

        xNeighbor = [ ...
            max(1,xx-1), ...
            min(imageWidth,xx+1)];

        yNeighbor = [ ...
            max(1,yy-1), ...
            min(imageHeight,yy+1)];

        starImage(yy,xNeighbor,:) = ...
            repmat( ...
            reshape(haloColor,1,1,3), ...
            1,2,1);

        starImage(yNeighbor,xx,:) = ...
            repmat( ...
            reshape(haloColor,1,1,3), ...
            2,1,1);
    end
end

end