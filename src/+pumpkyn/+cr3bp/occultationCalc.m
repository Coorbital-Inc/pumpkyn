function tauDur = occultationCalc(tau,x,rP1,rP2,muStar)
%% Purpose:
%
%  This routine will determine the occultations of Primary 1 as a result of
%  a satellite behind Pimary 2 in the CR3BP frame of reference.
%
%% Inputs:
%
%  tau                      [N x 1]                 Dimensionless Time
%
%  x                        [N x 6]                 Dimensionless Pos
%                                                   and Velocity of
%                                                   satellite in the CR3BP
%
%  rP1                      double                  Dimensionless Radius
%                                                   of Primary 1
%
%  rP2                      double                  Dimensionless Radius
%                                                   of Primary 2
%
%  muStar                   double                  Dimenionless Mass Ratio
%                                                   of both primaries
%                                                   muStar = m2/(m1+m2)
%
%% Outputs:
%
%  tauDur                   [M x 2]                 Intervals of
%                                                   dimensionless time
%                                                   where P2 is blocking
%                                                   or partially blocking
%                                                   the LOS of the
%                                                   satellite to P1
%% Revision History:
%  Darin C. Koblick                                         (c) 10-03-2025
%  Copyright 2025 Coorbital, Inc.
%% --------------------- Begin Code Sequence ------------------------------
if nargin == 0
                        tau0 = 2*pi;
                          Np = 15;
                          %% 
                          pm = +1;
[tau0, x0, mu, tStar, lStar] = pumpkyn.cr3bp.getTulip(tau0,Np,pm,1e-12);
                     [tau,x] = pumpkyn.cr3bp.prop(tau0,x0,mu);
                         rP2 = 1737.1./lStar;
                         rP1 = 6378.0./lStar;
                      tauDur = pumpkyn.cr3bp.occultationCalc(tau,x,rP1,rP2,mu);
                      occDur = diff(tauDur,1,2).*tStar./3600;  %Hrs
                      totOcc = sum(occDur);
                  return;
end

%% Get angle between P1 and P2 rel to SV:
     rS2P2 = [1-muStar,0,0] - x(:,1:3);  % SV->P2
     rS2P1 = [-muStar,0, 0] - x(:,1:3);  % SV->P1
     theta = pumpkyn.util.bsxAng(rS2P2,rS2P1,2);
     P2Ang = atand(rP2./pumpkyn.util.vmag(rS2P2,2)); 
     P1Ang = atand(rP1./pumpkyn.util.vmag(rS2P1,2));
    minAng = P1Ang + P2Ang;

%% Look for zero crossings:
             phi = theta - minAng;  
          phi_pm = zeros(size(phi));
 phi_pm(phi < 0) = -1;
 phi_pm(phi > 0) = +1;

%% Find the intervals corresponding to zero crossings:
      phi_pm_dot = diff(phi_pm,1,1);
          idxPos = find(phi_pm_dot > 0);
          idxNeg = find(phi_pm_dot < 0);
       tauPos_lb = tau(idxPos);
       tauPos_ub = tau(idxPos+1);
       tauNeg_lb = tau(idxNeg+1);
       tauNeg_ub = tau(idxNeg);
      tauPos_bnd = [tauPos_lb,tauPos_ub];
      tauNeg_bnd = [tauNeg_lb,tauNeg_ub];
      tauPos_bnd = [min(tauPos_bnd,[],2),max(tauPos_bnd,[],2)];
      tauNeg_bnd = [min(tauNeg_bnd,[],2),max(tauNeg_bnd,[],2)];
          tauDur = [0 0];

%% Refine crossings w/ spline interpolation:
if ~isempty(tauPos_bnd)
    tauPos = zeros(size(tauPos_bnd,1),1);
    tauNeg = zeros(size(tauNeg_bnd,1),1);
    for ti=1:size(tauPos_bnd,1)
            tauPos(ti,1) = fminbnd(@(x)abs(interp1(tau,phi,x,'spline')), ...
                                   tauPos_bnd(ti,1),tauPos_bnd(ti,2));
            tauNeg(ti,1) = fminbnd(@(x)abs(interp1(tau,phi,x,'spline')), ...
                                   tauNeg_bnd(ti,1),tauNeg_bnd(ti,2));
    end
    
    if tauPos(1) < tauNeg(1)
         tauDur = [tauPos,tauNeg];
    else
         tauDur = [tauNeg,tauPos];
    end
end

end