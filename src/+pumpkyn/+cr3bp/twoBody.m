function rvDot = twoBody(t,rv,mu,dim3)
%% Purpose:
%
%  This routine will propagate a satellite using the two-body equations
%  of motion.  All inputs and outputs are dimensionalized.
%
%% Inputs:
%
%  t                        [N x 1]             Dimensionalized Time
%
%  rv                       [N x 6]             Dimensionalized position
%                                               and velocity (km, and km/s)
%
%  mu                       double              Gravitaional Constant of
%                                               primary body
%
%  dim3                     integer             Singleton Dimension
%                                               Specifier (which dimension
%                                               are pos/vel located?)
%
%% Outputs:
%
%  rvDot                    [N x 6]             Derivatives of rv wrt time
%                                               (km, and km/s)
%
%% Revision History:
%  Darin C. Koblick                                         (c) 10-01-2025
%  Copyright 2025 Coorbital, Inc.
%% -------------------------- Begin Code Sequence -------------------------
if nargin == 0
   mu = 398600.4418;
   rv0 = [7000,0,0,0,0,sqrt(mu/7000)];
   opts = odeset('absTol',1e-10,'relTol',1e-10);
   [t,rv] = ode113(@pumpkyn.cr3bp.twoBody,[0 86400],rv0,opts,mu,1);
   figure('Color',[1 1 1]);
   plot3(rv(:,1),rv(:,2),rv(:,3));
   grid on;
   axis equal;
   return;
end
   [rv,fSeq] =  fDim(rv,dim3);
       rvDot =  NaN(size(rv));
rvDot(:,1:3) =  rv(:,4:6);
rvDot(:,4:6) = -mu.*rv(:,1:3)./vmag(rv(:,1:3),2).^3;
       rvDot = eDim(rvDot,fSeq);
end