function lla = toLLA(pos,muStar,lStar)
%% Purpose:
%
%  This routine will take dimensionless states in the rotating barycentric 
%  frame and convert them to dimensionalized 
%  observer latitude, longitude, and altitude with respect to the 
%  Moon selenographic frame of reference.
%
%  Note this is a pseudo-transformation aimed at placing points
%  appropriately between the lunar surface and the CR3BP.
%
%% Inputs:
%
%  pos                    [N x 3]           Dimensionless Position [x,y,z]
%                                           in the CR3BP (rotating
%                                           barycentric frame)
%
%  muStar               double              Mass ratio parameter
%                                           mu = m2/(m1+m2)
%
%  lStar                double              Characteristic Length
%                                           in same units as altitude
%
%% Outputs:
%
%  lla                  [N x 3]             Dimensionalized Latitude
%                                           Longitude, and Altitude in
%                                           [rad,rad,km]
%
%% Revision History:
%  Darin C. Koblick                                         (c) 10/10/2025
%  Copyright 2025 Coorbital, Inc.
%% -------------------------- Begin Code Sequence -------------------------
if nargin == 0
       [lat,lon] = meshgrid(linspace(-pi/2+pi/25,pi/2-pi/25,25), ...
                            linspace(-pi,+pi-2*pi/25,25));
          lla_in = [lat(:),lon(:),lon(:).*0+15];
          muStar = 0.012150585609624;        %Mass ratio
           lStar = 389703.264829278;
             pos = pumpkyn.cr3bp.fromLLA(lla_in,muStar,lStar);
             lla = pumpkyn.cr3bp.toLLA(pos,muStar,lStar);
             figure('Color',[1 1 1]);
             plot(pumpkyn.util.vmag(lla_in-lla,2));
             grid on;
             ylabel('Conversion Error');
       return;
end
%Planet Radius Constants:
     rMoon  = 1737.1;
%Convert from barycentric to lunarcentric:
        rPA = pos - [1-muStar,0,0];
        rPA = rPA.*lStar;
%Convert from cartesian to spherical:    
[lon,lat,r] = cart2sph(rPA(:,1),rPA(:,2),rPA(:,3));
        lla = [lat,lon,r-rMoon];
end