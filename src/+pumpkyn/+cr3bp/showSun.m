function [h,globe] = showSun(jd0,lStar,muStar,M,hIn)
%% Purpose:
%
%  Place or update the Sun in the Earth-Moon rotating barycentric frame.
%  Position and radius are scaled using the supplied scene characteristic
%  length. Passing a globe previously returned by showSun updates only its
%  transform, avoiding repeated mesh and texture creation.
%
%% Inputs:
%
%   jd0                 double              Epoch used to find the sun
%                                           in an ECI coordinate frame
%                                           to be converted to the CR3BP.
%
%  lStar                double              Characteristic Length (km)
%
%  muStar               double              Mass Ratio of Primaries
%                                           muStar = mu2/(mu1+mu2)
%
%  M                    double              Characteristic Mass of
%                                           both primaries M = m1 + m2 (kg)
%
%  hIn                  handle              Optional figure, axes, or Sun
%                                           globe returned by showSun
%
%% Outputs:
%
%  h                    handle              Handle to current axes
%
%  globe                handle              Handle to Globe graphics
%
%% Revision History:
%  Darin C. Koblick                                              07/27/26
%  Copyright 2026 Coorbital, Inc.
%% --------------------------- Begin Code Sequence ------------------------

if nargin == 0
    muStar = 0.012150585609624;
     lStar = 389703.264829278;
         M = 5.9736E24 + 7.35E22; %Characteristic mass (kg)
       jd0 = pumpkyn.util.juliandate('01/01/2000 00:00:00');
    pumpkyn.cr3bp.showSun(jd0,lStar,muStar,M);
    return;
end

validateattributes( ...
    jd0,{'numeric'},{'real','finite','scalar'}, ...
    'pumpkyn.cr3bp.showSun','jd0');

validateattributes( ...
    lStar,{'numeric'},{'real','finite','scalar','positive'}, ...
    'pumpkyn.cr3bp.showSun','lStar');

validateattributes( ...
    muStar,{'numeric'}, ...
    {'real','finite','scalar','>',0,'<',1}, ...
    'pumpkyn.cr3bp.showSun','muStar');

validateattributes( ...
    M,{'numeric'},{'real','finite','scalar','positive'}, ...
    'pumpkyn.cr3bp.showSun','M');

updateExistingGlobe = ...
    nargin >= 5 && isgraphics(hIn,'surface');

if updateExistingGlobe
    globe = hIn;
    h = ancestor(globe,'axes');

    validSunGlobe = ...
        isa( ...
        globe.Parent, ...
        'matlab.graphics.primitive.Transform') && ...
        strcmp(globe.Tag,'pumpkynSunGlobe') && ...
        strcmp(globe.Parent.Tag,'pumpkynSunTransform');

    if ~validSunGlobe

        error( ...
            'pumpkyn:showSun:InvalidGlobe', ...
            ['The supplied surface must be a globe previously ', ...
            'returned by pumpkyn.cr3bp.showSun.']);
    end

    sunTransform = globe.Parent;
elseif nargin < 5 || isempty(hIn)
    hFigure = figure('Color',[0 0 0]);
    h = axes('Parent',hFigure,'Color','k');
elseif isgraphics(hIn,'axes')
    h = hIn;
elseif isgraphics(hIn,'figure')
    h = gca(hIn);
else
    error( ...
        'pumpkyn:showSun:InvalidParent', ...
        'hIn must be a figure, axes, or Sun globe surface.');
end

%% Determine the Sun position in the rotating frame
[rSunJ2K,vSunJ2K] = pumpkyn.util.planetPosVel( ...
    jd0,'Earth','Sun');

[sunState,~,epochLStar] = pumpkyn.cr3bp.fromJ2K( ...
    jd0,[rSunJ2K,vSunJ2K],muStar,M,1,2);

% fromJ2K uses the instantaneous Earth-Moon distance. Convert the
% Earth-relative result to the characteristic length used by this scene.
         rE = [-muStar,0,0];
sunPosition = (sunState(1,1:3)-rE).*epochLStar./lStar + rE;

%% Create the reusable Sun globe when needed
if ~updateExistingGlobe
    originalNextPlot = h.NextPlot;
    restoreNextPlot = onCleanup( ...
        @() set(h,'NextPlot',originalNextPlot));

    h.NextPlot = 'add';
    % Keep the base mesh in physical kilometres. The transform applies
    % the scene scaling and translation without changing the texture data.
    [~,globe] = pumpkyn.util.sun3D( ...
        [0 0 0],true,1,h);
    sunTransform = hgtransform( ...
        'Parent',h, ...
        'Tag','pumpkynSunTransform');
    globe.Parent = sunTransform;
    globe.Tag = 'pumpkynSunGlobe';
    clear restoreNextPlot;
    set(h,'Color','k');
    axis(h,'equal','vis3d','off');
    grid(h,'off');
end

%% Scale and position the Sun
sunScale = eye(4);
sunScale(1:3,1:3) = eye(3)./lStar;
sunTransform.Matrix = ...
    makehgtform('translate',sunPosition)*sunScale;
end
