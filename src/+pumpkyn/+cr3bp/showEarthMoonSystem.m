function [h,earthGlobe,moonGlobe,sunGlobe] = ...
        showEarthMoonSystem( ...
        jd0,lStar,muStar,M,hIn,varargin)
%% showEarthMoonSystem
%
% Create or update an Earth-Moon-Sun scene in the Earth-Moon rotating
% barycentric frame. Body meshes and textures are created only once in an
% axes. Later calls with the same figure or axes update the Earth and Sun
% transforms at the requested Julian date while reusing the static Moon.
%
%% Inputs:
%
%  jd0                  double              Julian date (days)
%
%  lStar                double              Scene characteristic length
%                                           (km)
%
%  muStar               double              Primary mass ratio
%                                           mu2/(mu1+mu2)
%
%  M                    double              Combined primary mass (kg)
%                                           Optional; defaults to the
%                                           toolkit Earth and Moon masses
%
%  hIn                  handle              Optional figure or axes
%
%  Quality              char/string         'interactive' (default) or
%                                           'high'
%
%  ShowSun              logical             Show/update the Sun (default
%                                           true for the Julian-date API)
%
%% Outputs:
%
%  h                    handle              Scene axes
%
%  earthGlobe           surface             Earth surface
%
%  moonGlobe            surface             Moon surface
%
%  sunGlobe             surface             Sun surface, or empty when
%                                           ShowSun is false
%
%% Calling examples:
%
%  [ax,earth,moon,sun] = ...
%      pumpkyn.cr3bp.showEarthMoonSystem( ...
%      jd0,lStar,muStar,M);
%
%  % Reuse all existing graphics at a later epoch.
%  pumpkyn.cr3bp.showEarthMoonSystem( ...
%      jd1,lStar,muStar,M,ax);
%
%  The legacy call showEarthMoonSystem(lStar,muStar) remains supported and
%  preserves its original Earth-Moon-only behavior.
%
%% Revision History:
%  Darin C. Koblick                                         08/26/2025
%  Darin C. Koblick                                         07/27/2026
%  Copyright 2026 Coorbital, Inc.
%% -------------------------- Begin Code Sequence -------------------------

constants = pumpkyn.util.getConst('astro');
defaultMass = ...
    constants.Earth.Mass + constants.EarthMoon.Mass;
defaultJulianDate = pumpkyn.util.juliandate( ...
    '01/01/2000 00:00:00');

selfTest = nargin == 0;
legacySyntax = false;

if nargin == 0
    % Use a near-full-Moon geometry so the Sun is close to the
    % Earth-Moon line while remaining visible beside the Earth.
    jd0 = pumpkyn.util.juliandate('01/21/2000 06:40:00');
    lStar = 389703.264829278;
    muStar = 0.012150585609624;
    M = defaultMass;
elseif nargin == 1
    lStar = 389703.264829278;
    muStar = 0.012150585609624;
    M = defaultMass;
elseif nargin == 2
    % Backward compatibility: showEarthMoonSystem(lStar,muStar)
    legacySyntax = true;
    muStar = lStar;
    lStar = jd0;
    jd0 = defaultJulianDate;
    M = defaultMass;
elseif nargin == 3 && isgraphics(muStar)
    % Backward compatibility: showEarthMoonSystem(lStar,muStar,parent)
    legacySyntax = true;
    hIn = muStar;
    muStar = lStar;
    lStar = jd0;
    jd0 = defaultJulianDate;
    M = defaultMass;
elseif nargin == 4 && isgraphics(M)
    hIn = M;
    M = defaultMass;
elseif nargin < 3
    error( ...
        'pumpkyn:showEarthMoonSystem:NotEnoughInputs', ...
        ['Use showEarthMoonSystem(jd0,lStar,muStar) or the ', ...
        'legacy showEarthMoonSystem(lStar,muStar) syntax.']);
elseif nargin < 4 || isempty(M)
    M = defaultMass;
end

hasParent = exist('hIn','var') && ...
    ~isempty(hIn) && ...
    (isgraphics(hIn,'figure') || isgraphics(hIn,'axes'));

if exist('hIn','var') && ~isempty(hIn) && ~hasParent
    varargin = [{hIn},varargin];
end

validateattributes( ...
    jd0,{'numeric'},{'real','finite','scalar'}, ...
    'pumpkyn.cr3bp.showEarthMoonSystem','jd0');

validateattributes( ...
    lStar,{'numeric'},{'real','finite','scalar','positive'}, ...
    'pumpkyn.cr3bp.showEarthMoonSystem','lStar');

validateattributes( ...
    muStar,{'numeric'}, ...
    {'real','finite','scalar','>',0,'<',1}, ...
    'pumpkyn.cr3bp.showEarthMoonSystem','muStar');

validateattributes( ...
    M,{'numeric'},{'real','finite','scalar','positive'}, ...
    'pumpkyn.cr3bp.showEarthMoonSystem','M');

parser = inputParser;
parser.FunctionName = 'pumpkyn.cr3bp.showEarthMoonSystem';

addParameter( ...
    parser,'Quality','interactive', ...
    @(value) (ischar(value) && isrow(value)) || ...
    (isstring(value) && isscalar(value)));

addParameter( ...
    parser,'ShowSun',~legacySyntax, ...
    @(value) (islogical(value) && isscalar(value)) || ...
    (isnumeric(value) && isscalar(value) && ismember(value,[0 1])));

parse(parser,varargin{:});

quality = validatestring( ...
    parser.Results.Quality, ...
    {'interactive','high'}, ...
    parser.FunctionName, ...
    'Quality');

showSun = logical(parser.Results.ShowSun);

%% Resolve the target axes

if ~hasParent
    hFigure = figure('Color','k');
    h = axes('Parent',hFigure,'Color','k');
elseif isgraphics(hIn,'axes')
    h = hIn;
    hFigure = ancestor(h,'figure');
else
    hFigure = hIn;

    if isempty(hFigure.CurrentAxes)
        h = axes('Parent',hFigure,'Color','k');
    else
        h = hFigure.CurrentAxes;
    end
end

hFigure.CurrentAxes = h;

earthGlobe = findBody(h,'pumpkynEarthGlobe');
moonGlobe = findBody(h,'pumpkynMoonGlobe');
sunGlobe = findBody(h,'pumpkynSunGlobe');

initializeAxes = ...
    isempty(earthGlobe) && ...
    isempty(moonGlobe) && ...
    isempty(sunGlobe);

originalNextPlot = h.NextPlot;
restoreNextPlot = onCleanup( ...
    @() set(h,'NextPlot',originalNextPlot));

h.NextPlot = 'add';

%% Create or update Earth

if isempty(earthGlobe)
    [~,earthGlobe] = pumpkyn.cr3bp.showEarth( ...
        jd0,lStar,muStar,hFigure, ...
        'Quality',quality);

    earthGlobe.Tag = 'pumpkynEarthGlobe';
else
    pumpkyn.cr3bp.showEarth( ...
        jd0,lStar,muStar,earthGlobe, ...
        'Quality',quality);
end

%% Create the static Moon once

if isempty(moonGlobe)
    [~,moonGlobe] = pumpkyn.cr3bp.showMoon( ...
        lStar,muStar,h, ...
        'Quality',quality);

    moonGlobe.Tag = 'pumpkynMoonGlobe';
end

%% Create or update Sun

if showSun
    if isempty(sunGlobe)
        [~,sunGlobe] = pumpkyn.cr3bp.showSun( ...
            jd0,lStar,muStar,M,h);
    else
        pumpkyn.cr3bp.showSun( ...
            jd0,lStar,muStar,M,sunGlobe);

        sunGlobe.Visible = 'on';
    end
elseif ~isempty(sunGlobe)
    sunGlobe.Visible = 'off';
    sunGlobe = gobjects(0);
end

clear restoreNextPlot;

%% Configure only a newly initialized scene

set(h,'Color','k','Clipping','off');

if initializeAxes
    axis(h,'equal','vis3d','off');
    grid(h,'off');
end

if selfTest
    earthPosition = [-muStar,0,0];
    moonPosition = [1-muStar,0,0];

    % Look back toward the Earth from just behind the Moon. The small
    % cross-track offset prevents the Moon from fully occulting Earth.
    set(h,'Projection','perspective');
    campos(h,moonPosition + [0.4,0.012,0]);
    camtarget(h,earthPosition);
    camup(h,[0,0,1]);
    camva(h,5);

    % Keep the renderer's depth range local to the Earth-Moon system.
    % Clipping is disabled above, so the distant Sun remains visible
    % without forcing the nearby bodies out of the depth buffer.
    xlim(h,[earthPosition(1)-0.04,moonPosition(1)+0.06]);
    ylim(h,[-0.05,0.05]);
    zlim(h,[-0.05,0.05]);
    drawnow;
end

end


function globe = findBody(ax,tag)
%% findBody
% Return at most one valid body surface with the requested tag.

globe = findobj( ...
    ax, ...
    'Type','surface', ...
    'Tag',tag);

if isempty(globe)
    globe = gobjects(0);
else
    globe = globe(1);
end

end
