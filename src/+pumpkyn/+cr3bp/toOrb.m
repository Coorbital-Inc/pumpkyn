function oev = toOrb(tau,rvND,M,muStar,tStar,lStar,P)
%% Purpose:
%
%  This routine will take the dimensionless states in the rotating 
%  barycentric frame (CR3BP) and convert them to keplerian orbital 
%  elements with respect to either primary body (e.g. Earth or Moon) 
%
%% Inputs:
%
%  tau                  [N x 1]             Dimensionless Time
%
%
%  rvND                 [N x 6]             Dimensionless States in
%                                           the CR3BP
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
%
%  lStar                double              Characteristic Length [km]
%
%
%  P                    integer             Reference Primary Body 
%                                           1 = Earth Centered
%                                           2 = Moon Centered
%% Outputs:
%
%  oev                  [N x 6]             Keplerian Orbital Elements:
%                                           ---------
%                                           sma  (km)
%                                           ecc
%                                           inc  (rad)
%                                           argp (rad)
%                                           raan (rad)
%                                           trua (rad)
%
%% Revision History:
%  Darin C. Koblick                                         (c) 12/10/2025
%  Copyright 2025 Coorbital, Inc.
%% -------------------------- Begin Code Sequence -------------------------
   G = 6.67384e-20;                  % Gravitational constant (km^3/kg/s^2)
%Compute the Mass of the two primaries:
if P == 2
    mu = G * muStar * M;             % Standard Gravitaional Param of Moon
else
    mu = G * (1 - muStar) * M;       % Standard Gravitaional Param of Earth
end
%Convert from CR3BP to PCI: 
     rv = pumpkyn.cr3bp.toPCI(tau,rvND(:,1:6),muStar,tStar,lStar,P);
%Convert to Keplerian Orbital Elements:
    oev = pumpkyn.cr3bp.eci2orb(mu,rv(:,1:3),rv(:,4:6),2);
end