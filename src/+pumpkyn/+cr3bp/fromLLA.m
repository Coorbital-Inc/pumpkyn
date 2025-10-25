function pos = fromLLA(lla,muStar,lStar)
%% Purpose:
%
%  This routine will take observer latitude, longitude, and altitude
%  with respect to the Moon (selenographic) and convert to dimensionless 
%  states in the rotating barycentric frame such that they can be used 
%  for calculations within the CR3BP.
%
%  Note this is not the same thing as going from a Principal Axis to
%  J2000 to CR3BP. It is a pseudo-transformation aimed at placing points
%  appropriately in the CR3BP.
%
%
%% Inputs:
%
%  lla                  [N x 3]             Dimensionalized Latitude
%                                           Longitude, and Altitude in
%                                           [rad,rad,km]
%
%  muStar               double              Mass ratio parameter
%                                           mu = m2/(m1+m2)
%
%  lStar                double              Characteristic Length
%                                           in same units as altitude
%
%% Outputs:
%
%  pos                    [N x 3]           Dimensionless Position [x,y,z]
%                                           in the CR3BP (rotating
%                                           barycentric frame)
%
%% Revision History:
%  Darin C. Koblick                                         (c) 10/10/2025
%  Copyright 2025 Coorbital, Inc.
%% -------------------------- Begin Code Sequence -------------------------
if nargin == 0
       [lat,lon] = meshgrid(linspace(-pi/2+pi/25,pi/2-pi/25,25), ...
                            linspace(-pi,+pi-2*pi/25,25));
             lla = [lat(:),lon(:),lon(:).*0+15];
          muStar = 0.012150585609624;        %Mass ratio
           lStar = 389703.264829278;
             pos = pumpkyn.cr3bp.fromLLA(lla,muStar,lStar);
       pumpkyn.cr3bp.showMoon()
       plot3(pos(:,1),pos(:,2),pos(:,3),'.r','MarkerSize',15);
       axis equal;
       return;
end
%Planet Radius Constants:
 rMoon  = 1737.1;
[x,y,z] = sph2cart(lla(:,2),lla(:,1),rMoon+lla(:,3));
    rPA = [x,y,z]./lStar;
    pos = rPA + [1-muStar,0,0];
end