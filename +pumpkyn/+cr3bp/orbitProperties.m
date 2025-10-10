function data = orbitProperties(x0,tau0,muStar,lStar)
%% Purpose:
%
%  This routine efficiently computes and aggregates key diagnostic 
%  quantities for a periodic orbit in the Circular Restricted 
%  Three-Body Problem (CR3BP). Given an initial state, period, and mass 
%  ratio, the function evaluates key parameters such as the Jacobi constant, 
%  pseudo-potential, stability index, and other invariants. The results 
%  are returned in a structured format to support both quick diagnostics 
%  and mission design trade studies.
%
%
%% Inputs:
%
%  x0                   [1 x 6]             Dimensionless States (pos/vel)
%
%  tau0                 [1 x 1]             Dimensionless Period
%
%  muStar               double              Mass Ratio of Primaries
%                                           muStar = mu2/(mu1+mu2)
%
%% Outputs:
%
%  data                 struct              Structure containing key
%                                           orbit properties:
%                                           -----------------
%                                           - Jacobi (Constant)           x
%                                           - StabilityIndex              x
%                                           - PseudoPotential             x
%                                           - Perilune (Min Dist From P2) x
%                                           - Apolune  (Max Dist From P2) x
%                                           - Perigee  (Min Dist from P1) x
%                                           - Apogee   (Max Dist from P1) x
%                                           - DoublingTime (ND)           x
%                                           - CharExponent                x
%                                           - StableEigenVals             x
%                                           - UnstableEigenVals           x
%                                           - StableEigenVecs             x
%                                           - UnstableEigenVecs           x
%                                           - MaxLunarOcc
%                                           - TotLunarOcc
%                                                                                         
%% Revision History:
%  Darin C. Koblick                                         (c) 10/2/2025
%  Copyright 2025 Coorbital, Inc.
%% ------------------------- Begin Code Sequence --------------------------
if nargin == 0
       x0 = [8.2339081983651485E-1	-1.9017764504099543E-28	9.8941366235910004E-4	...
            -2.3545391932685812E-15	1.2634272983881797E-1	2.2367029429442455E-16];
     tau0 = 2.7430007981241529E+0;
   muStar = 1.215058560962404E-2;
    lStar = 389703;
     data = pumpkyn.cr3bp.orbitProperties(x0,tau0,muStar,lStar);
       %1.1804065333857600E+3
       return;
end

%% Propagate Oribt a full Period
           [tau,x] = pumpkyn.cr3bp.prop(tau0,[x0(:); reshape(eye(6),[36 1])],muStar);

%% Extract the Monodromy Matrix:
                 M = reshape(x(end,7:end),[6 6]);

%% Find the Eigen Values and Vectors:
             [V,D] = eig(M);
            lambda = diag(D);

%% Compute the Stability Index:
data.StabilityIndex = max(abs(0.5.*(lambda + 1./lambda)));

%% Compute the Jacobi Constant and Pseudo Potential:
[data.Jacobi,data.PseudoPotential] = pumpkyn.cr3bp.jacobi(x0,muStar);

%% Find Perilune and Apolune (ND) when dr/dt = 0:
[data.Perilune,data.Apolune] = minMax(x,[1-muStar,0,0,0,0,0]);

%% Find Perigee and Apogee (ND) when dr/dt = 0:
  [data.Perigee,data.Apogee] = minMax(x,[-muStar,0,0,0,0,0]);
    
%% Find Stable and Unstable Eigen Values/Vectors:
                       idx_u = abs(lambda) > (1 + 1e-3);
      data.UnstableEigenVals = lambda(idx_u);   %[N x 1]
      data.UnstableEigenVecs = V(:,idx_u)';     %[N x 6]
        data.StableEigenVals = lambda(~idx_u);  %[M x 1]
        data.StableEigenVecs = V(:,~idx_u)';    %[M x 6]

%% Doubling Time and Charateristic Exponent:
           data.CharExponent = NaN;
           data.DoublingTime = NaN;
     if data.StabilityIndex >= 1.01
        data.CharExponent = max(real((1./tau0).*log(data.UnstableEigenVals))); 
        data.DoublingTime = min(real(tau0.*log(2)./log(data.UnstableEigenVals)));
     end

%% Compute Lunar Occultations
if exist('lStar','var')
              rMoon = 1737.1./lStar;
             rEarth = 6378.0./lStar;
             tauOcc = pumpkyn.cr3bp.occultationCalc(tau,x,rEarth,rMoon,muStar);
             tauDur = diff(tauOcc,1,2);
   data.MaxLunarOcc = max(tauDur);
   data.TotLunarOcc = sum(tauDur);
end

end

function [rMin,rMax,rMag] = minMax(rv,rvBody)
%% Purpose:
%
%  Find the global min and global max values
%
                  rvSat2B = rv(:,1:6) - rvBody;
                     rMag = pumpkyn.util.vmag(rvSat2B(:,1:3),2);
                     rMin = min(rMag); 
                     rMax = max(rMag);
               
end
