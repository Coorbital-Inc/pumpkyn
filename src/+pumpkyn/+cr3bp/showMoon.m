function [h,globe,hIn] = showMoon(lStar,muStar,hIn,varargin)
%% Purpose:
%
%  This routine will properly place the moon in dimensionless coordinates
%  with the correct scaling and position at [1-mu,0,0]
%
%% Inputs:
%
%  lStar                double              Characteristic Length (km)
%
%  muStar               double              Mass Ratio of Primaries
%                                           muStar = mu2/(mu1+mu2)
%
%  hIn                  handle              Optional handle input
%
%  Quality              char                'high' (default) or
%                                           'interactive'
%
%% Outputs:
%
%  h                    handle              Handle to current axes
%
%  globe                handle              Handle to Globe graphics
%
%% Revision History:
%  Darin C. Koblick                                              08/26/2025
%  Copyright 2025 Coorbital, Inc.
%% --------------------------- Begin Code Sequence ------------------------

if nargin == 0
        lStar = 389703.264829278;
       muStar = 0.012150585609624;
    [h,globe] = pumpkyn.cr3bp.showMoon(lStar,muStar);
    return;
end
hasParent = nargin >= 3 && ...
    ~(ischar(hIn) || (isstring(hIn) && isscalar(hIn)));

if nargin >= 3 && ~hasParent
    varargin = [{hIn},varargin];
end

quality = 'interactive';

if ~isempty(varargin)
    parser = inputParser;
    parser.FunctionName = 'pumpkyn.cr3bp.showMoon';

    addParameter( ...
        parser,'Quality',quality, ...
        @(value) (ischar(value) && isrow(value)) || ...
        (isstring(value) && isscalar(value)));

    parse(parser,varargin{:});

    quality = validatestring( ...
        parser.Results.Quality, ...
        {'interactive','high'}, ...
        parser.FunctionName, ...
        'Quality');
end

if ~hasParent
    hIn = figure('color',[0 0 0]);
    set(gca(hIn),'color','k');
end

if isgraphics(hIn,'axes')
    ax = hIn;
else
    ax = gca(hIn);
end

[h,globe] = pumpkyn.util.moon3D( ...
    [1-muStar,0,0],true,1/lStar,ax, ...
    'Quality',quality);
set(ax,'color','k','clipping','off');
axis(ax,'off','equal');
end
