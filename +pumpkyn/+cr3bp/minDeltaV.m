function dV = minDeltaV(x0,x1,muStar)
%% Purpose:
%
%  This routine will compute the theoretical minimum cost associated with
%  changing the energy of the spacecraft from an initial orbit to a final
%  orbit without considering its transfer path.
%
%% Inputs:
%
%  x0                   [N x 6]             Dimensionless States (pos/vel)
%                                           of first orbit
%
%  x1                   [N x 6]             Dimensionless States (pos/vel)
%                                           of second orbit
%
%  muStar               double              Mass Ratio of Primaries
%                                           muStar = mu2/(mu1+mu2)
%
%% Outputs:
%
% dV                    [N x 1]             Minimum Velocity required
%                                           to change energy from x0 to x1
%                                          (dimensionless)
%
%% Revision History:
%  Darin C. Koblick                                              09-24-2025
%  Copyright 2025 Coorbital, Inc.
%% --------------------------- Begin Code Sequence ------------------------
if nargin == 0
                        Np = 7;               
                      tau0 = (5/6)*2*pi;
                        pm = -1;                      
[~,r,v,muStar,lStar,tStar] = pumpkyn.cr3bp.tulipConstellation(Np,tau0,1,pm,0);
                    x0 = [r(1,:),v(1,:)];
                  tau1 = 2*pi;
               [~,r,v] = pumpkyn.cr3bp.tulipConstellation(Np,tau1,1,pm,0);
                    x1 = [r(1,:),v(1,:)];
                    dV = pumpkyn.cr3bp.minDeltaV(x0,x1,muStar);
                dV_kms = dV.*lStar./tStar;
                disp(dV_kms);
                    return;
end
%Compute the jacobi constant and pseudo-potential for x0:
[J0,U0] = pumpkyn.cr3bp.jacobi(x0,muStar);
[J1,U1] = pumpkyn.cr3bp.jacobi(x1,muStar);
%Min dV is found when U* is maximized:
  Ustar = max([U0,U1],[],2);
Ustar_2 = 2.*Ustar;
     dV = abs(sqrt(Ustar_2 - J1) - sqrt(Ustar_2 - J0));
end