function [dV_LB,dV_UB] = stationKeeping_deltaV(x0,tau0,dVMaxErr,dTau,muStar,N)
%% Purpose:
%
%  This routine will compute the upper bound associated with the station
%  keeping delta-v costs required to maintain a periodic orbit. 
%
%% Inputs:
%
%  x0                   [1 x 6]             Dimensionless States (pos/vel)
%
%  tau0                 double              Dimensionless Period of Orbit
%
%  dVMaxErr             double              Maximum dispersion velocity
%                                           error (dimensionless)
%
%  dTau                 double              Time between correction
%                                           maneuvers
%                                           (dimensionless)
%
%  muStar               double              Mass ratio of primaries
%                                           mu = m2/(m1+m2)
%
%  N                    double              Number of burns/year
%
%% Outputs:
%
%  dV_LB                double              Total annual Delta-V lower
%                                           bound (dimensionlesS)
%
%  dV_UB                double              Total annual Delta-V upper
%                                           bound (dimensionless)
%
%% Revision History:
%  Darin C. Koblick                                         (c) 09/30/2025
%  Copyright 2025 Coorbital, Inc.
%% ------------------------ Begin Code Sequence ---------------------------
if nargin == 0
                          Np = 13; %Number of petals
                        tau0 = 2*pi; 
                          pm = +1;
[tau0,x0,muStar,tStar,lStar] = pumpkyn.cr3bp.getTulip(tau0,Np,pm,1e-12);
                dVMaxErr = sqrt(3)*(0.02/1000)*tStar/lStar; % 1cm/s -> ND
                    dTau = (1:1:60)'.*86400/tStar; %Correction interval
                       N = (365.25*86400/tStar)./dTau;
                    dV_LB = NaN(size(N));
                    dV_UB = NaN(size(N));
          for tt=1:numel(dTau)            
            [dV_LB(tt,1),dV_UB(tt,1)]  = pumpkyn.cr3bp.stationKeeping_deltaV(x0,tau0, ...
                                         dVMaxErr,dTau(tt),muStar,N(tt));
          end
                   dV_LB = (dV_LB.*lStar./tStar).*1000; %m/s
                   dV_UB = (dV_UB.*lStar./tStar).*1000; %m/s    
          figure('Color',[1 1 1]);
          plot(dTau.*tStar./86400,dV_LB,'b'); hold on;
          plot(dTau.*tStar./86400,dV_UB,'r');
          grid on;
          xlabel('\Deltat_b [days]'); ylabel('\DeltaV_{tot} [m/s/yr]');

return;
end
         ubm = 3.0;          %Factor of 3 to account for nav and exec error
         
        data = pumpkyn.cr3bp.orbitProperties(x0,tau0,muStar);

        if isempty(data.UnstableEigenVals)
            %No instabilities
            [dV_LB,dV_UB] = deal(0);
            return;
        end

%% Compute Velocity Error after dTau:
   dV_b_LB = dVMaxErr.*2.^(dTau./data.DoublingTime);

%% Determine Number of burn corrections per year:
   dV_b_UB = ubm.*dV_b_LB;   
     dV_UB = dV_b_UB.*N;   %Total dV per calendar year  (Upper bound)
     dV_LB = dV_b_LB.*N;   %Total dV per calendar year  (lower bound)
end