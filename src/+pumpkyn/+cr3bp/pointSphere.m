function [r,lla,xyz] = pointSphere(N,dr,rP)
%% Purpose:
%
%  The purpose of this routine is to create a sphere of points that could
%  resemble observers on the second primary body. This point sphere can be
%  used for discrete GDOP and LOS analysis.
%
%% Inputs:
%
%  N                Integer                     Number of Unique 
% %                                             Points in sphere
%
%  
%  dr               [1 x 3]                     Dimensionless Offset of
%                                               sphere e.g.  [1-mu,0,0]
%                                               for second primary body
%   
%  rP               double                      Dimensionless Radius of
%                                               the sphere
%
%
%% Outputs:
%
%  r               [M x 3]                      Dimensionless position of
%                                               each point on the sphere
%
% lla              [M x 3]                      lat,lon,alt of each point
%                                               on sphere, 0 deg longitude
%                                               corresponds to -x
%                                               direction (rad,rad,ND)
%
% xyz              [P x Q x 3]                  Dimensionless position
%                                               of each point on the
%                                               sphere in a format
%                                               compatable with surf
%
%% Revision History:
%  Darin C. Koblick                                              07/30/2025
%  Copyright 2025 Coorbital, Inc.
%% ------------------------ Begin Code Sequence ---------------------------

if nargin == 0
   
        N = 1500;
       mu = 1.215058560962404E-2;
    lStar = 389703.264829278;
       rP = 1738.1./lStar;
       dr = [1-mu,0,0];
  [r,lla] = pumpkyn.cr3bp.pointSphere(N,dr,rP);
      idx = (abs(lla(:,2)) < 10*pi/180 | abs(lla(:,2)) > 350*pi/180) & abs(lla(:,1)) < 10*pi/180;
      
    figure('color',[1 1 1]);
    plot3(r(:,1),r(:,2),r(:,3),'.k'); hold on;
    plot3(r(idx,1),r(idx,2),r(idx,3),'.r');
    %plot3(-mu,0,0,'.b','markersize',20);
    axis equal;
    return;
end
     [~,~,~,xS,yS,zS] = pumpkyn.util.sphere3D(N);
                 x = xS(:);
                 y = yS(:);
                 z = zS(:);
           [az,el] = cart2sph(x,y,z);
                az = mod(az-pi,2*pi);
                 r = [x,y,z].*rP;
                 r = r + dr;
               lla = [el,az,zeros(size(el))];
               xyz = cat(3,reshape(r(:,1),size(xS)), ...
                           reshape(r(:,2),size(xS)), ...
                           reshape(r(:,3),size(xS)));
end
