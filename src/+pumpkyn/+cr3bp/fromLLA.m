function x = fromLLA(jd,lla,muStar,lStar,P)
%% Purpose:
%
%  This routine will take observer latitude, longitude, and altitude
%  with respect to either the Earth (geodetic) or the Moon (selenographic) 
%  and convert to dimensionless states in the rotating barycentric 
%  frame such that they can be used for calculations within the CR3BP.
%
%  Note this is not the same thing as going from a Principal Axis to
%  J2000 to CR3BP. It is a pseudo-transformation aimed at placing points
%  appropriately in the CR3BP.
%
%
%% Inputs:
%
%  jd                   [N x 1]             Julian Date in Days
%
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
%  P                    integer             Primary Body Specifier:
%                                           1 = Earth
%                                           2 = Moon
%
%% Outputs:
%
%  x                     [N x 6]           Dimensionless Position [x,y,z]
%                                          in the CR3BP (rotating
%                                          barycentric frame) and velocity
%                                          in the CR3BP frame
%
%% Revision History:
%  Darin C. Koblick                                         (c) 06/18/2026
%  Copyright 2026 Coorbital, Inc.
%% -------------------------- Begin Code Sequence -------------------------
if nargin == 0
       [lat,lon] = meshgrid(linspace(-pi/2+pi/25,pi/2-pi/25,25), ...
                            linspace(-pi,+pi-2*pi/25,25));
             lla = [lat(:),lon(:),lon(:).*0+15];
          muStar = 0.012150585609624;        %Mass ratio
           lStar = 389703.264829278;
              jd = [];
             pos = pumpkyn.cr3bp.fromLLA(jd,lla,muStar,lStar,2);
       pumpkyn.cr3bp.showMoon()
       plot3(pos(:,1),pos(:,2),pos(:,3),'.r','MarkerSize',15);
       

       %Verify specific marker on near side (Sea of Rains):
       pos = pumpkyn.cr3bp.fromLLA(jd,[32.8*pi/180 -15.6*pi/180 0],muStar,lStar,2);
       plot3(pos(:,1),pos(:,2),pos(:,3),'.g','MarkerSize',15);
       %Verify specific marker on far side (Mare Moscoviense):
       pos = pumpkyn.cr3bp.fromLLA(jd,[27.3*pi/180  147.9*pi/180 0],muStar,lStar,2);
       plot3(pos(:,1),pos(:,2),pos(:,3),'.b','MarkerSize',15);
       axis equal;

       return;
end

if P == 2
         rMoon  = 1737.1;
        [x,y,z] = sph2cart(lla(:,2),lla(:,1),rMoon+lla(:,3));
            rPA = [x,y,z]./lStar;
            pos = rPA + [1-muStar,0,0];
            vel = pos.*0;
              x = [pos,vel];
else
              M = 5.9736E24 + 7.35E22;              %Characteristic mass
    [rECI,vECI] = pumpkyn.util.LLA2ECI(jd,lla,2);
              x = pumpkyn.cr3bp.fromJ2K(jd,[rECI,vECI],muStar,M,1);
end

end