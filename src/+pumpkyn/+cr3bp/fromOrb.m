function rv = fromOrb(t,oev,M,muStar,tStar,lStar,P)
%% Purpose:
%
%  This routine will take keplerian orbital elements with respect to
%  either primary body (e.g. Earth or Moon) and convert them to 
%  dimensionless states in the rotating barycentric frame such that they 
%  can be used with the CR3B equations of motion.
%
%% Inputs:
%
%  t                    [N  x 1]            P = 1: Julian date [days]
%                                           P = 2: Time from epoch [s]
%
%  oev                  [N x 6]             Keplerian Orbital Elements:
%                                           ---------
%                                           sma  (km)
%                                           ecc
%                                           inc  (rad)
%                                           argp (rad)
%                                           raan (rad)
%                                           trua (rad)
%                                           P = 1: J2000 equatorial axes
%                                           P = 2: Ideal CR3BP inertial axes
%
%  M                    double              Characteristic mass of
%                                           primaries M = m1 + m2
%                                           in kg
%
%  muStar               double              Mass ratio parameter
%                                           mu = m2/(m1+m2)
%
%  tStar                double              Characteristic Time [s]
%
%  lStar                double              Characteristic Length [km]
%
%  P                    integer             Reference Primary Body 
%                                           1 = Earth Centered
%                                           2 = Moon Centered
%% Outputs:
%
%  rv                   [N x 6]             Dimensionless States in
%                                           the CR3BP
%
%% Revision History:
%  Darin C. Koblick                                         (c) 12/10/2025
%  Copyright 2025 Coorbital, Inc.
%% -------------------------- Begin Code Sequence -------------------------

   G = 6.67384e-20;                % Gravitational constant (km^3/kg/s^2)
  
%Compute the Mass of the two primaries:
if P == 2
    mu = G * muStar * M;             % Standard Gravitaional Param of Moon
else
    mu = G * (1 - muStar) * M;       % Standard Gravitaional Param of Earth
end

%Convert from orbital elements to PCI:
[r0,v0] = pumpkyn.cr3bp.orb2eci(mu, oev, 2);

%Convert to the appropriate rotating reference frame:
if P == 1
    [rv,tStar2,lStar2] = pumpkyn.cr3bp.fromJ2K(t(:),[r0,v0],muStar,M,P,2);
              rEarth = [-muStar,0,0];
          rv(:,1:3) = (rv(:,1:3)-rEarth).*lStar2./lStar + rEarth;
          rv(:,4:6) = rv(:,4:6).*(lStar2./tStar2).*(tStar./lStar);
else
     rv = pumpkyn.cr3bp.fromPCI(t./tStar,[r0,v0],muStar,tStar,lStar,P);
end
end
