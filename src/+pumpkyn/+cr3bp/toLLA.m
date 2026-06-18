function lla = toLLA(jd,pos,muStar,lStar,P)
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
%  jd                   [N x 1]             Julian Date in Days
%
%
%  pos                  [N x 3]             Dimensionless Position [x,y,z]
%                                           in the CR3BP (rotating
%                                           barycentric frame)
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
               P = 1;
              jd = pumpkyn.util.juliandate('03/21/2030 00:00:00'); 
             pos = pumpkyn.cr3bp.fromLLA(jd,lla_in,muStar,lStar,P);
             lla = pumpkyn.cr3bp.toLLA(jd,pos,muStar,lStar,P);
             err = lla_in-lla;
             idx = err(:,2) > pi;
             err(idx,2) = err(idx,2)-2*pi;
             idx = err(:,2) < -pi;
             err(idx,2) = err(idx,2)+2*pi;
             figure('Color',[1 1 1]);
             plot(pumpkyn.util.vmag(err,2));
             grid on;
             ylabel('Conversion Error');
       return;
end

if P == 2
    %Planet Radius Constants:
         rMoon  = 1737.1;
    %Convert from barycentric to lunarcentric:
            rPA = pos(:,1:3) - [1-muStar,0,0];
            rPA = rPA.*lStar;
    %Convert from cartesian to spherical:    
    [lon,lat,r] = cart2sph(rPA(:,1),rPA(:,2),rPA(:,3));
            lla = [lat,lon,r-rMoon];
else
              M = 5.9736E24 + 7.35E22;              %Characteristic mass
              x = pumpkyn.cr3bp.toJ2K(jd,pos,muStar,M,1);
            lla = pumpkyn.util.ECI2LLA(jd,x(:,1:3),2);
end

end