function [v0,vf,converged] = minLambert(rv0,rvf,tof,muStar,tStar,lStar,M,dim3)
%% Purpose:
%
%  This routine will compute a single-revolution direct lambert-style
%  transfer betwen two points in the CR3BP. A two-body Lambert solution 
%  (via pyKep) seeds an initial guess, and convergence is achieved using 
%  the linearized dynamics from the propagated State Transition Matrix in
%  the form of a targeting algorithm.
%
%  Instead of needing the user to guess the transfer direction through
%  use of a prograde/retrograde specifier, this routine will compute the
%  transfer both ways and output the solution associated with the minimum
%  total delta-V. For stability purposes, it will also consider transfers
%  using both the larger and smaller primary for computing the two-body
%  lambert transfer.
%
%% Inputs:
%
%  rv0                      [1 x 6]             Initial Pos/Vel at t=0
%                                               (Dimensionless)
%
%  rvf                      [1 x 6]             Final Pos/Vel at t=tof
%                                               (Dimensionless)
%
%  tof                      double              Dimensionless Time of
%                                               Flight
%
%  muStar                   double              Mass Ratio of Primaries
%                                               muStar = m2/(m1+m2)
%
%  tStar                    double              Characteristic Time (s)
%
%
%  lStar                    double              Characteristic Length (km)
%
%
%
%  M                        double              Characteristic Mass of
%                                               both Primaries (kg)
%                                               M = m1 + m2
%
%  dim3                     integer             Singleton Dimension
%                                               Specifier
%
%% Outputs:
%
%  v0                       [1 x 3]             Initial Dimensionless 
%                                               Velocity at r0
%                                               corresponding to t = 0
%
%  vf                       [1 x 3]             Final Dimensinless Velocity
%                                               at rf corresponding to
%                                               t = tof
%
%% Revision History:
%  Darin C. Koblick                                         (c) 10/23/2025
%  Copyright 2025 Coorbital, Inc.
%% -------------------------- Begin Code Sequence -------------------------

if nargin == 0
    muStar = 0.012150585609624;
     lStar = 389703.264829278;
     tStar = 382981.289129055;
         M = 5.9736E24 + 7.35E22;  %Characteristic mass (sum of primaries)
       rv0 = [0.98841303973042, 0, 0.03950808689588, 0, 0, 0];
       rvf = [0.955615424345515,-0.0155940896976095, -0.0120526582538886, 0, 0, 0];
       tof = 0.5*86400/tStar; %10*86400/tStar;
       [v0,vf,converged] = pumpkyn.cr3bp.minLambert(rv0,rvf,tof, ...
                           muStar,tStar,lStar,M,2);
       %Show trajectory:
             rv0(:,4:6) = v0;
          [~,rv] = pumpkyn.cr3bp.prop(tof,rv0,muStar);
        figure('Color',[1 1 1]);
        plot3(rv(:,1),rv(:,2),rv(:,3)); hold on;
        plot3(rv0(1,1),rv0(1,2),rv0(1,3),'.g','markersize',12);
        plot3(rvf(1,1),rvf(1,2),rvf(1,3),'.r','markersize',12);
        axis equal;
        grid on;
       return;
end

%% Flatten input arrays:
[rv0,fSeq] = pumpkyn.util.fDim(rv0,dim3);
       rvf = pumpkyn.util.fDim(rvf,dim3);

%% Compute Lambert Transfer about both primaries (P1 and P2):
[converged_P1,v0_P1,vf_P1,dVtot_P1] = minLambertSol(rv0,rvf,tof, ...
                                                    muStar,tStar,lStar,M,1);
[converged_P2,v0_P2,vf_P2,dVtot_P2] = minLambertSol(rv0,rvf,tof, ...
                                                    muStar,tStar,lStar,M,2);

%% Determine Minimum dVtot Solution:
if converged_P1 && (~converged_P2 || dVtot_P1 < dVtot_P2)
converged = converged_P1;
       v0 = v0_P1;
       vf = vf_P1;
elseif converged_P2
converged = converged_P2;
       v0 = v0_P2;
       vf = vf_P2;
else
converged = false;
       v0 = rv0(:,4:6);
       vf = rvf(:,4:6);
end

%% Reshape outputs:
fSeq.postShift(2) = 3;
v0 = pumpkyn.util.eDim(v0,fSeq);
vf = pumpkyn.util.eDim(vf,fSeq);
end

function [converged,v0,vf,dVtot] = minLambertSol(rv0,rvf,tof,muStar,tStar,lStar,M,P)
%% Compute Lambert Transfer both ways:
[v0_p,vf_p,converged_p] = pumpkyn.cr3bp.directLambert(rv0(:,1:3), ...
                          rvf(:,1:3),tof,+1,muStar,tStar,lStar,M,P,2);

[v0_m,vf_m,converged_m] = pumpkyn.cr3bp.directLambert(rv0(:,1:3), ...
                          rvf(:,1:3),tof,-1,muStar,tStar,lStar,M,P,2);


%% Determine the total dV associated w/ each direction:
dV1_p = v0_p - rv0(:,4:6);
dV2_p = rvf(:,4:6) - vf_p;
 dV_p = norm(dV1_p) + norm(dV2_p);

dV1_m = v0_m - rv0(:,4:6);
dV2_m = rvf(:,4:6) - vf_m;
 dV_m = norm(dV1_m) + norm(dV2_m);

%% Select the output w/ the lowest overall dV: 
   v0 = v0_p;
   vf = vf_p;
dVtot = dV_p;

if dV_m < dV_p
     v0 = v0_m;
     vf = vf_m;
  dVtot = dV_m;
end

%% Total Convergence:
converged = converged_p & converged_m;

end