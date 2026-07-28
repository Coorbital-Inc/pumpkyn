function [r,v] = planetPosVel(jd0,center,target,varargin)
%% Purpose:
%
%  This routine computes the position and velocity of one celestial body
%  with respect to another. It can use the high-fidelity planetEphemeris
%  model or the faster analytical Sun and Moon models included with
%  Pumpkyn. Dense planetEphemeris requests are evaluated on a reduced
%  Julian Date grid and interpolated back to the requested epochs.
%
%% Calling Syntax:
%
%  [r,v] = pumpkyn.util.planetPosVel(jd0,center,target)
%
%  [r,v] = pumpkyn.util.planetPosVel( ...
%      jd0,center,target,'Method',method)
%
%% Inputs:
%
%  jd0                  [N x 1]             Julian Date Vector (Days)
%
%  center               char/string         Reference body or point from
%                                           which the position and velocity
%                                           of the target body are computed
%
%                                           Options include:
%                                           'Sun', 'Mercury', 'Venus',
%                                           'Earth', 'Moon', 'Mars',
%                                           'Jupiter', 'Saturn', 'Uranus',
%                                           'Neptune', 'Pluto',
%                                           'SolarSystem', and 'EarthMoon'
%
%  target               char/string         Target body or point of interest
%                                           whose position and velocity are
%                                           computed with respect to center
%
%                                           Options include:
%                                           'Sun', 'Mercury', 'Venus',
%                                           'Earth', 'Moon', 'Mars',
%                                           'Jupiter', 'Saturn', 'Uranus',
%                                           'Neptune', 'Pluto',
%                                           'SolarSystem', and 'EarthMoon'
%
%  Method               char/string         Ephemeris calculation method:
%
%                                            'auto' (default) uses
%                                            planetEphemeris when available
%                                            and otherwise uses an analytical
%                                            model for supported body pairs.
%
%                                            'toolbox' requires and uses
%                                            planetEphemeris.
%
%                                            'analytic' uses the faster
%                                            Pumpkyn Sun or Moon model. This
%                                            method supports Earth-Moon and
%                                            Earth-Sun states in either
%                                            direction.
%
%% Outputs:
%
%  r                    [N x 3]             Position of target with respect
%                                           to center
%                                           [km,km,km]
%
%  v                    [N x 3]             Velocity of target with respect
%                                           to center
%                                           [km/s,km/s,km/s]
%
%% Revision History:
%  Darin C. Koblick                                         (c) 06/24/2026
%
%  Copyright 2026 Coorbital, Inc.
%% -------------------------- Begin Code Sequence -------------------------
if nargin == 0
    center = 'Earth';
    target = 'Moon';
    jd0 = pumpkyn.util.juliandate('01/01/2000 12:00:00');
      t = (0:24*60:86400*365)';
     jd = jd0 + t./86400;
    tTime = tic;
    [r,v] = pumpkyn.util.planetPosVel( ...
        jd,center,target,'Method','auto');
    eTime = toc(tTime);
    tTime = tic;
[r_i,v_i] = planetEphemeris(jd,center,target);
   eiTime = toc(tTime);
       dr = pumpkyn.util.vmag(r-r_i,2);
       dv = pumpkyn.util.vmag(v-v_i,2);

    figure('color',[1 1 1]);
    subplot(1,2,1);
    plot(t./86400,100.*dr./pumpkyn.util.vmag(r_i,2));
    grid on;
    ylabel('Position Error (RSS) [%]');
    xlabel('Time [Days]');
    subplot(1,2,2);
    plot(t./86400,100.*dv./pumpkyn.util.vmag(v_i,2));
    grid on;
    ylabel('Velocity Error (RSS) [%]');
    xlabel('Time [Days]');

    fprintf('moonPosVel self-test:\n');
    fprintf('  Execution Speed-up: %.6g \n',eiTime./eTime);
    fprintf('  Max position error: %.6g km\n',max(dr));
    fprintf('  RMS position error: %.6g km\n',rms(dr));
    fprintf('  Max velocity error: %.6g m/s\n',max(dv).*1000);
    fprintf('  RMS velocity error: %.6g m/s\n',rms(dv).*1000);
    [r,v] = deal([]);
    return;
end

method = parseMethod(varargin);

% Tuning Parameters:
dtGridDays   = 6/24;                    % Discrete steps
minInterpPts = 4;                       % Min number of points to interp
padDays      = minInterpPts*dtGridDays; % Padding for interpolation
     iMethod = 'spline';                % Interpolation Method
% Input Checks:
if isempty(jd0)
    r = zeros(0,3);
    v = zeros(0,3);
    return;
end
jd0 = jd0(:);

useAnalytic = strcmp(method,'analytic');

if ~useAnalytic
    hasPlanetEphemeris = ...
        exist('planetEphemeris','file') ~= 0;

    if ~hasPlanetEphemeris
        if strcmp(method,'auto')
            useAnalytic = true;
        else
            error( ...
                'pumpkyn:planetPosVel:PlanetEphemerisUnavailable', ...
                ['Method=''toolbox'' requires planetEphemeris, but ', ...
                 'planetEphemeris is not available on the MATLAB path.']);
        end
    end
end

if useAnalytic
    [r,v] = evaluateAnalyticModel(jd0,center,target);
    return;
end

% Min/Max of requested timespan
jdMin = min(jd0);
jdMax = max(jd0);

if isscalar(jd0) || jdMax == jdMin
    [r1,v1] = planetEphemeris(jd0(1),center,target);
          r = repmat(r1,numel(jd0),1);
          v = repmat(v1,numel(jd0),1);
else
    jdStart = jdMin - padDays;
    jdStop  = jdMax + padDays;
      nGrid = max(minInterpPts,ceil((jdStop - jdStart)/dtGridDays) + 1);
      jd0_i = linspace(jdStart,jdStop,nGrid).';
  [r_i,v_i] = planetEphemeris(jd0_i,center,target);
          r = interp1(jd0_i,r_i,jd0,iMethod);
          v = interp1(jd0_i,v_i,jd0,iMethod);
end

end

function method = parseMethod(optionArguments)
%% Parse the Method name-value option.

method = 'auto';

if isempty(optionArguments)
    return;
end

if numel(optionArguments) ~= 2 || ...
        ~isTextScalar(optionArguments{1}) || ...
        ~strcmpi(char(optionArguments{1}),'Method')

    error( ...
        'pumpkyn:planetPosVel:InvalidOptions', ...
        ['Specify the optional calculation method as ', ...
         '''Method'',''auto'', ''Method'',''toolbox'', or ', ...
         '''Method'',''analytic''.']);
end

methodValue = optionArguments{2};

if ~isTextScalar(methodValue)
    error( ...
        'pumpkyn:planetPosVel:InvalidMethod', ...
        'Method must be ''auto'', ''toolbox'', or ''analytic''.');
end

method = validatestring( ...
    char(methodValue), ...
    {'auto','toolbox','analytic'}, ...
    'pumpkyn.util.planetPosVel', ...
    'Method');

end


function [r,v] = evaluateAnalyticModel(jd0,center,target)
%% Evaluate a supported analytical relative-state model.

if strcmpi(center,'Earth') && strcmpi(target,'Moon')
    [r,v] = pumpkyn.util.moonPosVel(jd0);
elseif strcmpi(center,'Moon') && strcmpi(target,'Earth')
    [r,v] = pumpkyn.util.moonPosVel(jd0);
    r = -r;
    v = -v;
elseif strcmpi(center,'Earth') && strcmpi(target,'Sun')
    [r,v] = pumpkyn.util.sunPosVel(jd0);
elseif strcmpi(center,'Sun') && strcmpi(target,'Earth')
    [r,v] = pumpkyn.util.sunPosVel(jd0);
    r = -r;
    v = -v;
else
    error( ...
        'pumpkyn:planetPosVel:FallbackUnsupported', ...
        ['Method=''analytic'' supports only Earth-Moon and ', ...
         'Earth-Sun states in either direction.']);
end

end


function valid = isTextScalar(value)
%% Return true for a character row or scalar string.

valid = ...
    (ischar(value) && isrow(value)) || ...
    (isstring(value) && isscalar(value));

end
