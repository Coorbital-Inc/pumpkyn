function [v0,vf,converged] = directLambert(r0,rf,tof,pm,muStar,tStar,lStar,M,P,dim3)
%% Purpose:
%
%  This routine will compute a single-revolution direct lambert-style
%  transfer betwen two points in the CR3BP. A two-body Lambert solution 
%  (via pyKep) seeds an initial guess, and convergence is achieved using 
%  the linearized dynamics from the propagated State Transition Matrix in
%  the form of a targeting algorithm.
%
%  For multiple revolution transfers, consider using global non-linear
%  optimization or breaking up the transfer into multiple points and
%  optimizing the results of this directLambert routine.
%
%  Important Limitations: 
%  - The STM is a local linear model valid only for small perturbations 
%    around the reference trajectory.
%
%  - If your transfer spans regions where gravity from each primary 
%    changes substantially (e.g., near Lagrange points), the linear 
%    approximation can break down quickly.
%
%  - In CR3BP, trajectories diverge exponentially near unstable 
%    manifolds—so a linear model's predictive power fades rapidly.
%
%% Inputs:
%
%  r0                       [1 x 3]             Initial Position at t=0
%                                               (Dimensionless)
%
%  rf                       [1 x 3]             Final Position at t=tof
%                                               (Dimensionless)
%
%  tof                      double              Dimensionless Time of
%                                               Flight
%
%  pm                       integer             +1 = prograde
%                                               -1 = retrograde
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
%  P                        integer             Which Primary has
%                                               the strongest gravitational
%                                               acceleration on SV 
%                                               throughout trajectory?
%                                               P = 1, (Earth)
%                                               P = 2, (Moon)
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
%  Darin C. Koblick                                         (c) 10/22/2025
%  Copyright 2025 Coorbital, Inc.
%% -------------------------- Begin Code Sequence -------------------------
if nargin == 0
    muStar = 0.012150585609624;
     lStar = 389703.264829278;
     tStar = 382981.289129055;
         M = 5.9736E24 + 7.35E22;  %Characteristic mass (sum of primaries)
        r0 = [0.98841303973042, 0, 0.03950808689588];
        rf = [0.955615424345515,-0.0155940896976095, -0.0120526582538886];
        pm = +1;
       tof = 0.5*86400/tStar; %10*86400/tStar;
         P = 2;
       [v0,vf,converged] = pumpkyn.cr3bp.directLambert(r0,rf,tof,pm, ...
                           muStar,tStar,lStar,M,P,2);

       %Show trajectory:
             rv0 = [r0,v0];
          [~,rv] = pumpkyn.cr3bp.prop(tof,rv0,muStar);
        figure('Color',[1 1 1]);
        plot3(rv(:,1),rv(:,2),rv(:,3)); hold on;
        plot3(r0(1,1),r0(1,2),r0(1,3),'.g','markersize',12);
        plot3(rf(1,1),rf(1,2),rf(1,3),'.r','markersize',12);
        axis equal;
        grid on;
       
       return;
end

[r0,fSeq] = pumpkyn.util.fDim(r0,dim3);
       rf = pumpkyn.util.fDim(rf,dim3);

maxIter = 25;

%% Compute Standard Gravitational Parameters of Both Primaries (kg):
     G = 6.67384e-20;             
   mu1 = G.*(1-muStar).*M;
   mu2 = G.*muStar.*M;
    mu = mu2; 
  if P == 1
      mu = mu1;
  end

%% Convert from CR3BP to PCI @ tof:
   rv0 = pumpkyn.cr3bp.toPCI(0,[r0,r0.*0],muStar,tStar,lStar,P);
   rvf = pumpkyn.cr3bp.toPCI(tof,[rf,rf.*0],muStar,tStar,lStar,P);
 
%% Solve Initial Problem Using 2-body Lambert:                 
    v0 = pumpkyn.pykep.lambert2Body(rv0(:,1:3), rvf(:,1:3), ...
                                     tof*tStar, mu, pm, 0, 2);

%% Convert States back to CR3BP:
           rv0(:,4:6) = v0(1,:);
                rvND0 = pumpkyn.cr3bp.fromPCI(0,rv0,muStar,tStar,lStar,P);

%% Converge From Initial Lambert Solution to CR3BP Solution
                 PHI0 = eye(6);        
                 opts = odeset('RelTol',1e-10,'AbsTol',1e-10);
            converged = false;

            for i = 1:maxIter
            [~,rvPHI] = ode113(@pumpkyn.cr3bp.eom,[0 tof], ...
                        [rvND0,PHI0(:)'],opts,muStar,1);       
                 PHIf = reshape(rvPHI(end,7:42),[6 6]);
                rStar = rvPHI(end,1:3);
                   dR = (rf-rStar)';
                   dV = (PHIf(1:3, 4:6) \ dR)';
                    if norm(dR) < 1e-6 || norm(dV) < 1e-12
                        converged = true;
                        break;
                    end
                   %Update State for next iteration: 
                   %  disp(norm(dR))
                  rvND0(:,4:6) = rvND0(:,4:6) + dV;
                 %% 
            end

%% Reshape initial and final velocities:            
v0 = pumpkyn.util.eDim(rvND0(1,4:6),fSeq);
vf = pumpkyn.util.eDim(rvPHI(end,4:6),fSeq);
end