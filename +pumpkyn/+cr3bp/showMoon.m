function [h,globe] = showMoon(lStar,muStar,hIn)
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
    pumpkyn.cr3bp.showMoon(lStar,muStar);
    return;
end
if ~exist('hIn','var')
    hIn = figure('color',[0 0 0]);
    set(gca(hIn),'color','k');
end

[h,globe] = pumpkyn.util.moon3D([1-muStar,0,0],true,1/lStar,gca(hIn));
axis(h, 'off','equal');

end