function [r,v] = planetPosVel(jd0,center,target)
%% Purpose:
%
%  This routine will compute the position and velocity of one planet with
%  respect to another using planetEphemeris default ephemeris
%  model. The ephemeris is evaluated on a reduced Julian Date grid and
%  interpolated back to the requested epochs for improved execution speed.
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
      t = (0:60:86400*27)';
     jd = jd0 + t./86400;
    tTime = tic;
    [r,v] = pumpkyn.util.planetPosVel(jd,center,target);
    eTime = toc(tTime);
    tTime = tic;
[r_i,v_i] = planetEphemeris(jd,center,target);
   eiTime = toc(tTime);
       dr = pumpkyn.util.vmag(r-r_i,2);
       dv = pumpkyn.util.vmag(v-v_i,2);

    figure('color',[1 1 1]);
    subplot(1,2,1);
    plot(t./86400,dr);
    grid on;
    ylabel('Position Error (RSS) [km]');
    xlabel('Time [Days]');
    subplot(1,2,2);
    plot(t./86400,dv);
    grid on;
    ylabel('Velocity Error (RSS) [km/s]');
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

hasPlanetEphemeris = exist('planetEphemeris','file') ~= 0;

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
% Min/Max of timespan
  jd0 = jd0(:);
jdMin = min(jd0);
jdMax = max(jd0);

%Call vectorized fallbacks:
if ~hasPlanetEphemeris
   if strcmpi(center,'Earth') && strcmpi(target,'Moon')
                    [r,v] = pumpkyn.util.moonPosVel(jd0);
   elseif strcmpi(center,'Moon') && strcmpi(target,'Earth')
                    [r,v] = pumpkyn.util.moonPosVel(jd0);
                        r = -r;
                        v = -v;
   else
   error('pumpkyn:planetPosVel:FallbackUnsupported', ...
         ['planetEphemeris is unavailable. The moonPosVel ', ...
          'fallback supports only Earth-to-Moon and ', ...
          'Moon-to-Earth states.']);
   end
   return;
end

% Input Checks:
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