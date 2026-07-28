function [hLogoAxes,hLogo] = addFigureLogo( ...
    fig,logoFile,logoPosition)
%% Purpose:
%
% Add a fixed-size PNG logo to the upper-left corner of a figure.
%
% logoPosition is the initial normalized [left bottom width height]
% bounding box. The logo is fitted within this box without changing its
% aspect ratio. Its size and upper-left margins are converted to pixels and
% remain fixed when the figure is resized.
%

%% Self test

if nargin == 0
    fig = figure( ...
        'Color','k', ...
        'Position',[100 100 800 600], ...
        'Name','addFigureLogo Self Test', ...
        'NumberTitle','off');

    logoFile = 'logo.png';
    logoPosition = [0.01 0.94 0.18 0.05];

    [hLogoAxes,hLogo] = pumpkyn.util.addFigureLogo( ...
        fig,logoFile,logoPosition);

    return;
end

%% Defaults and validation

if nargin < 2 || isempty(logoFile)
    logoFile = 'logo.png';
end

if nargin < 3 || isempty(logoPosition)
    logoPosition = [0.01 0.94 0.18 0.05];
end

if ~isgraphics(fig,'figure')
    error('fig must be a valid figure handle.');
end

if ~isnumeric(logoPosition) || ...
        ~isequal(size(logoPosition),[1 4]) || ...
        any(logoPosition(3:4) <= 0)

    error( ...
        'logoPosition must be [left bottom width height].');
end

%% Resolve the image path

logoFile = char(logoFile);
utilityFolder = fileparts(mfilename('fullpath'));

% A bare default filename always means the packaged Pumpkyn wordmark.
% This prevents an unrelated logo.png in the current folder from being
% selected accidentally.
[logoFolder,logoName,logoExtension] = fileparts(logoFile);
isDefaultLogo = ...
    isempty(logoFolder) && ...
    strcmpi([logoName logoExtension],'logo.png');

if isDefaultLogo
    logoFile = fullfile(utilityFolder,'logo.png');
elseif ~isfile(logoFile)
    logoFile = fullfile(utilityFolder,logoFile);
end

if ~isfile(logoFile)
    error('Logo file not found: %s',logoFile);
end

%% Preserve the active axes and replace an existing logo

activeAxes = fig.CurrentAxes;

oldLogoAxes = findall( ...
    fig,'Type','axes','Tag','figureLogoAxes');

if ~isempty(oldLogoAxes)
    delete(oldLogoAxes);
end

%% Read the PNG

[logoImage,~,logoAlpha] = imread(logoFile);

if ~isempty(logoAlpha) && isinteger(logoAlpha)
    logoAlpha = double(logoAlpha)./ ...
        double(intmax(class(logoAlpha)));
end

% Remove transparent padding so the requested size applies to the visible
% wordmark and MATLAB does not downsample unused pixels.
if ~isempty(logoAlpha)
    [visibleRows,visibleColumns] = find(logoAlpha > 0);

    if ~isempty(visibleRows)
        rowRange = min(visibleRows):max(visibleRows);
        columnRange = min(visibleColumns):max(visibleColumns);

        logoImage = logoImage(rowRange,columnRange,:);
        logoAlpha = logoAlpha(rowRange,columnRange);
    end
end

%% Convert the initial normalized placement to fixed pixels

figPixels = getpixelposition(fig);
figWidth = figPixels(3);
figHeight = figPixels(4);

maximumSize = max(1,round([ ...
    logoPosition(3)*figWidth, ...
    logoPosition(4)*figHeight]));

imageSize = size(logoImage);
imageAspectRatio = imageSize(2)/imageSize(1);

if maximumSize(1)/maximumSize(2) > imageAspectRatio
    logoHeight = maximumSize(2);
    logoWidth = round(logoHeight*imageAspectRatio);
else
    logoWidth = maximumSize(1);
    logoHeight = round(logoWidth/imageAspectRatio);
end

% geometry = [leftMargin topMargin width height]
geometry = round([ ...
    logoPosition(1)*figWidth, ...
    (1-logoPosition(2)-logoPosition(4))*figHeight, ...
    max(1,logoWidth), ...
    max(1,logoHeight)]);

logoPixelPosition = [ ...
    geometry(1), ...
    figHeight-geometry(2)-geometry(4), ...
    geometry(3), ...
    geometry(4)];

%% Create and draw the logo

hLogoAxes = axes( ...
    'Parent',fig, ...
    'Units','pixels', ...
    'Position',logoPixelPosition, ...
    'Color','none', ...
    'Visible','off', ...
    'HitTest','off', ...
    'PickableParts','none', ...
    'HandleVisibility','off', ...
    'Clipping','off', ...
    'Tag','figureLogoAxes');

hLogo = image( ...
    hLogoAxes,logoImage, ...
    'HitTest','off', ...
    'PickableParts','none', ...
    'Interpolation','bilinear');

if ~isempty(logoAlpha)
    hLogo.AlphaData = logoAlpha;
    hLogo.AlphaDataMapping = 'none';
end

axis(hLogoAxes,'image');
axis(hLogoAxes,'off');

% image can reset its parent axes when NextPlot is 'replace'. Apply the
% overlay-only behavior after drawing so the logo never regains a toolbar
% or competes with the foreground scene for mouse input.
set(hLogoAxes, ...
    'Color','none', ...
    'Visible','off', ...
    'HitTest','off', ...
    'PickableParts','none', ...
    'HandleVisibility','off', ...
    'Clipping','off', ...
    'Tag','figureLogoAxes');

if isprop(hLogoAxes,'Toolbar')
    hLogoAxes.Toolbar = [];
end

if isprop(hLogoAxes,'Interactions')
    hLogoAxes.Interactions = [];
end

uistack(hLogoAxes,'top');

%% Install a figure-specific resize callback

previousResizeFcn = fig.SizeChangedFcn;
previousAutoResize = fig.AutoResizeChildren;

fig.AutoResizeChildren = 'off';

fig.SizeChangedFcn = @(src,event) resizeLogo( ...
    src,event,previousResizeFcn,hLogoAxes,geometry);

% Restore the original figure behavior when this logo is deleted.
hLogoAxes.DeleteFcn = @(~,~) restoreResizeBehavior( ...
    fig,previousResizeFcn,previousAutoResize);

%% Restore the previously active axes

if ~isempty(activeAxes) && isgraphics(activeAxes,'axes')
    fig.CurrentAxes = activeAxes;
end

drawnow;

end

function resizeLogo(fig,event,previousFcn,hLogoAxes,geometry)
%% Run the existing resize callback, then reposition the fixed-size logo.

runCallback(previousFcn,fig,event);

if ~isgraphics(hLogoAxes,'axes')
    return;
end

figPixels = getpixelposition(fig);

hLogoAxes.Position = [ ...
    geometry(1), ...
    figPixels(4)-geometry(2)-geometry(4), ...
    geometry(3), ...
    geometry(4)];

end

function runCallback(callback,source,event)
%% Execute a previously installed MATLAB callback.

if isempty(callback)
    return;
elseif isa(callback,'function_handle')
    callback(source,event);
elseif iscell(callback)
    feval(callback{1},source,event,callback{2:end});
else
    evalin('base',char(callback));
end

end

function restoreResizeBehavior(fig,resizeFcn,autoResize)
%% Restore the figure configuration when the logo is deleted.

if isgraphics(fig,'figure') && ...
        strcmp(fig.BeingDeleted,'off')

    fig.SizeChangedFcn = resizeFcn;
    fig.AutoResizeChildren = autoResize;
end

end
