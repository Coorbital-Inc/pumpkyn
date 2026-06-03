function [h,globe] = showSun(psi0,lStar,hIn)
%% Purpose:
%
%  This routine will properly place the Sun in dimensionless coordinates
%  with the correct scaling and position defined by alpha0.
%
%% Inputs:
%
%   psi0                double              Sun angle in Radians
%                                           with respect to the E-M system
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
%  Darin C. Koblick                                              05/20/2026
%  Copyright 2026 Coorbital, Inc.
%% --------------------------- Begin Code Sequence ------------------------

if nargin == 0
      psi0 = 0;
     lStar = 389703.264829278;
    pumpkyn.cr3bp.showSun(psi0,lStar);
    return;
end
if ~exist('hIn','var')
    hIn = figure('color',[0 0 0]);
    set(gca(hIn),'color','k');
end
     %Sun-to-barycenter vector:
     rSun = 389.1723985.*[cos(psi0), sin(psi0), 0];
[h,globe] = pumpkyn.util.sun3D(rSun,true,1/lStar,gca(hIn));

end