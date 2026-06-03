function [h,globe,moon] = showEarthMoonSystem(lStar,muStar,h)
%% Purpose:
%
%  This routine will properly place the moon in dimensionless coordinates
%  with the correct scaling and position at [1-mu,0,0] and the Earth
%  in dimensionless coordinates with correct scaling and position at
%  [-mu,0,0].
%
%% Inputs:
%
%  lStar                double              Characteristic Length (km)
%
%  muStar               double              Mass Ratio of Primaries
%                                           muStar = mu2/(mu1+mu2)
%
%   h                   handle              Optional figure handle input
%
%% Outputs:
%
%  h                    handle              Handle to current axes
%
%  globe                handle              Handle to Globe graphics
%
%  moon                 handle              Handle to Moon graphics
%
%% Revision History:
%  Darin C. Koblick                                              08/26/2025
%  Copyright 2025 Coorbital, Inc.
%% --------------------------- Begin Code Sequence ------------------------

if nargin == 0
     lStar = 389703.264829278;
    muStar = 0.012150585609624;
    [h,globe,moon] = pumpkyn.cr3bp.showEarthMoonSystem(lStar,muStar);
    return;
end

if ~exist('hIn','var')
      h = figure('color',[0 0 0]);
    set(gca(h),'color','k');
end

[h,globe] = pumpkyn.cr3bp.showEarth(lStar,muStar,h);
 [h,moon] = pumpkyn.cr3bp.showMoon(lStar,muStar,h);

 set(gca,'Clipping','off');
end