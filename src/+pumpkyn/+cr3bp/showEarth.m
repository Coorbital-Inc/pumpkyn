function [h,globe] = showEarth(jd0,lStar,muStar,hIn)
%% Purpose:
%
%  This routine will properly place the Earth in dimensionless coordinates
%  with the correct scaling and position at [-mu,0,0]
%
%% Inputs:
%
%  jd0                  double              Julian Date (days)
%
%
%  lStar                double              Characteristic Length (km)
%
%  muStar               double              Mass Ratio of Primaries
%                                           muStar = mu2/(mu1+mu2)
%
%  hIn                  handle              Optional handle input
%
%% Outputs:
%
%  h                    handle              Handle to current axes
%
%  globe                handle              Handle to Globe graphics
%
%% Revision History:
%  Darin C. Koblick                                              06/18/2026
%  Copyright 2026 Coorbital, Inc.
%% --------------------------- Begin Code Sequence ------------------------
if nargin == 0
          jd0 = pumpkyn.util.juliandate('03/21/2000 00:00:00');
        lStar = 389703.264829278;
       muStar = 0.012150585609624;
    [h,globe] = pumpkyn.cr3bp.showEarth(jd0,lStar,muStar);
    return;
end

if ~exist('hIn','var')
    hIn = figure('color',[0 0 0]);
    set(gca(hIn),'color','k');
    axis(gca(hIn), 'off','equal');
    hold on;
end
             opts.posOffset = [-muStar,0,0];
                 opts.scale = 1./lStar;  
                  opts.type = 'day';
               opts.overlay = true;
                 opts.atmos = false; 
            opts.AddShading = false;
                 opts.WGS84 = true;

[h,globe] = pumpkyn.util.earth3D(opts,hIn);

if ~isempty(jd0)
    % Get three Earth-fixed basis directions expressed in CR3BP frame
    xObs_x = pumpkyn.cr3bp.fromLLA(jd0,[0 0 0],muStar,lStar,1);
    xObs_y = pumpkyn.cr3bp.fromLLA(jd0,[0 pi/2 0],muStar,lStar,1);
    xObs_z = pumpkyn.cr3bp.fromLLA(jd0,[pi/2  0 0],muStar,lStar,1);
    % Convert absolute positions into Earth-centered direction vectors
    rE = [-muStar, 0, 0];
    ex = xObs_x(1:3) - rE;
    ey = xObs_y(1:3) - rE;
    ez = xObs_z(1:3) - rE;
    ex = ex ./ norm(ex);
    ey = ey ./ norm(ey);
    ez = ez ./ norm(ez);
    DCME2C = [ex(:), ey(:), ez(:)];
    xyzData = [globe.XData(:), globe.YData(:), globe.ZData(:)];
    xyzLocal = xyzData - rE;
    % Rotate the surface points on the globe:
    xyzData = pumpkyn.util.multiplyDCM(repmat(DCME2C,[1 1 size(xyzLocal,1)]), ...
        xyzLocal,2) + rE;
    %Update the surface points:
    globe.XData = reshape(xyzData(:,1),size(globe.XData));
    globe.YData = reshape(xyzData(:,2),size(globe.YData));
    globe.ZData = reshape(xyzData(:,3),size(globe.ZData));
end

set(gca(h),'clipping','off');
end