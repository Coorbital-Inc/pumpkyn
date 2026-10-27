function [yDot,Ht,dHdy,aThrust] = tfMinEoM(tau,y,Tmax,c,muStar) %#ok<INUSD>
%% Purpose:
%
%  This routine contains the low-thrust equations of motion in the CR3BP
%  designed to be integrated by a numerical integration routine such as
%  ODE78, ODE113, and ODE45.  This routine will also compute the
%  corresponding state transition matrix along with the states (optional).
%
%  References:
%  2015. AIAA. JGCD. Low-Thrust Minimum-Fuel Optimization in the Circular
%  Restricted Three-Body Problem.
%
%% Inputs:
%
%  tau                      double                  Dimensionless Time
%
%
%  y                        [14 x 1]                Initial States
%                                                   Dimensionless:
%                                                   -------
%                                                   r [x,y,z]
%                                                   v [dr/dt]
%                                                   m [ND]
%                                                   lambda_r [x,y,z]
%                                                   lambda_v [x,y,z]
%                                                   lambda_m 
%
%  Tmax                     double                  Maximum Thrust
%                                                   Magnitude
%
%  c                        double                  Exhaust Velocity
%                                                   (c = Isp*g0) of
%                                                   specific impulse
%
%  muStar                   double                  Mass Ratio of
%                                                   two primaries
%                                                   muStar = m2/(m1+m2)
%
%% Outputs:
%
%  yDot                     [14 x 1]                Derivative of Initial
%                                                   Dimensionless States
%                                                   With Respect to
%                                                   dimensionless time tau.
%% Revision History:
%  Darin C. Koblick                                          (c) 10/27/2025
%  Factored gravity / Hamiltonian derivatives                    10/19/2026
%  Copyright Coorbital Inc.
%% ------------------------ Begin Code Sequence ---------------------------
if nargin == 0
   %Table Three: Minimum Time Solutions Computed:
   muStar = 1.21506683E-2;
    lStar = 3.84405000E5;   %characteristic length (km)
    tStar = 3.75676967E5;   %characteristic time (s)
       m0 = 1500;           %Initial Satellite mass (kg)
       g0 = 9.80665*tStar*tStar/lStar;              %m/s^2 -> ND
     Tmax = 10;                                     %Newtons
     Tmax = (Tmax/m0)*tStar^2/(lStar*1000);         %Max Thrust ND 
     %  tf = 7.8549*86400/tStar;                     %Days ->ND
      Isp = 3000/tStar;    %Specific Impulse s -> ND
        c = Isp*g0;         %Exhuast Velocity (ND)

       %Boundary Conditions:
       r0 = [-0.019488511458668,-0.016033479812051,0];
       v0 = [8.918881923678198,-4.081793688818725,0];

       %Initial States:
       lambda_r0 = [15.616017, 32.875896, -0.094522];
       lambda_v0 = [-0.101606, 0.044791, -0.000150];
       lambda_m0 = 0.133266;
            PHI0 = eye(14);
        y = [r0(:); v0(:); 1; lambda_r0(:); lambda_v0(:); lambda_m0; PHI0(:)];
     yDot = pumpkyn.cr3bp.tfMinEoM(0,y,Tmax,c,muStar);
        return;
end
%Extract initial states:
        r = y(1:3);
        v = y(4:6);
        m = y(7);
 lambda_r = y(8:10);
 lambda_v = y(11:13);
 lambda_m = y(14);
lambda_vx = lambda_v(1);
lambda_vy = lambda_v(2);
lambda_vz = lambda_v(3);
     PHI = [];
if numel(y) > 14
     PHI = reshape(y(15:210),[14 14]);
end
%Re-assign variables:
       x = r(1);
       y = r(2);
       z = r(3);
    xDot = v(1); 
    yDot = v(2);
%Definition of g(r) and h(v) functions (Equation 2)
       d = sqrt((x + muStar)^2 + y^2 + z^2);
       r = sqrt((x - 1 + muStar)^2 + y^2 + z^2);
      d3 = d^3;
      r3 = r^3;
      d5 = d^5;
      r5 = r^5;
      d7 = d^7;
      r7 = r^7;
      gr = [x - (1-muStar)*(x+muStar)/d3 - muStar*(x-1+muStar)/r3
            y - (1-muStar)*y/d3 - muStar*y/r3
            0 - (1-muStar)*z/d3 - muStar*z/r3];
      hv = [+2.*yDot
            -2.*xDot;
            +0];
%See CR3BP Equations of Motion in eom.m:
%  d = sqrt((x + muStar)^2 + y^2 + z^2);
%  r = sqrt((x - 1 + muStar)^2 + y^2 + z^2);
% d3 = d^3;
% r3 = d^3;
%   xDDot =  x - (1-muStar)*(x+muStar)/d3 - muStar*(x-1+muStar)/r3 + 2*ydot;
%   yDDot =  y - (1-muStar)*y/d3 - muStar*y/r3 - 2*xdot;
%   zDDot =  -(1-muStar)*z/d3 - muStar*z/r3;

%Determine G = d g(r) / dr
% G = jacobian(gr,[x,y,z]);, see Eq 8
% clear d3 d5 r3 r5;
% G = subs(G,((muStar + x)^2 + y^2 + z^2)^(3/2),'d3');
% G = subs(G,((muStar + x)^2 + y^2 + z^2)^(5/2),'d5');
% G = subs(G,((muStar + x - 1)^2 + y^2 + z^2)^(3/2),'r3');
% G = subs(G,((muStar + x - 1)^2 + y^2 + z^2)^(5/2),'r5');
     G = zeros(3,3);
G(1,1) = (muStar - 1)/d3 - muStar/r3 - (3*(muStar + x)*(muStar + x)*(muStar - 1))/d5 + ...
         (3*muStar*(muStar + x - 1)*(muStar + x - 1))/r5 + 1;
G(1,2) = (3*muStar*y*(muStar + x - 1))/r5 - (3*y*(muStar + x)*(muStar - 1))/d5;
G(1,3) = (3*muStar*z*(muStar + x - 1))/r5 - (3*z*(muStar + x)*(muStar - 1))/d5;
G(2,1) = G(1,2);
G(2,2) = (muStar - 1)/d3 - muStar/r3 + (3*muStar*y^2)/r5 - (3*y^2*(muStar - 1))/d5 + 1;
G(2,3) = (3*muStar*y*z)/r5 - (3*y*z*(muStar - 1))/d5;
G(3,1) = G(1,3);
G(3,2) = G(2,3);
G(3,3) = (muStar - 1)/d3 - muStar/r3 + (3*muStar*z^2)/r5 - (3*z^2*(muStar - 1))/d5;
%Determine H = d h(v) / dv
%  H = jacobian(hv,[xDot,yDot,zDot]);, see Eq 8
         H = zeros(3,3);
    H(1,2) = 2;
    H(2,1) = -H(1,2);
%Determine the minimum-time switching function (Equation 19):
           S = -norm(lambda_v)*c/m - lambda_m;
           u = 1;
           if S > 0
                u = 0;
           end
%Hamiltonian for Minimum Time (Equation 17):
       alpha = -lambda_v./norm(lambda_v);                  %Equation 10
     aThrust =  u.*alpha.*Tmax./m;                         %Thrust Aceel Term
        rDot =  v;
        vDot =  gr + hv + aThrust;
        mDot = -u*Tmax/c;

if nargout > 1      
          Ht = 1 + dot(lambda_r,rDot,1) + dot(lambda_v,vDot,1) + lambda_m*mDot;
end

if nargout > 2
%Compute dHdy = [dH/dx, dH/dy, dH/dz, dH/dv, dH/dm, dH/dLr, dH/dLv, dHdLm]
    dHdy = [lambda_v'*G,lambda_r'+lambda_v'*H, ...
            Tmax*u*norm(lambda_v)/m^2,v',vDot',mDot];
end

%Costate Derivatives:
 lambdaDot_r = -G'*lambda_v;
 lambdaDot_v = -lambda_r - H'*lambda_v;
 lambdaDot_m = -norm(lambda_v)*u*Tmax/m^2;

%All State Derivatives (Equation 15):
        yDot = [rDot(:); 
                vDot(:); 
                mDot; 
                lambdaDot_r(:); 
                lambdaDot_v(:); 
                lambdaDot_m];

%Compute the Derivative of the State Transition Matrix (STM):
% jacobian(y, [r,v,m,lambda_r,lambda_v,lambda_m])
% A = [14 x 14]
if ~isempty(PHI)
            A = zeros(14,14);
 lambda_v_Mag = norm(lambda_v);
        A(1,:) = [0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0];
        A(2,:) = [0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0];
        A(3,:) = [0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0];
    A(4:6,1:3) =  G;
    A(4:6,4:6) =  H;
 A(8:10,11:13) = -G';
A(11:13,11:13) = -H';

  A(4,7) = Tmax*lambda_vx*u/(lambda_v_Mag*m^2);

 A(4,11) = -(Tmax*u*(lambda_vy^2 + lambda_vz^2))/(lambda_v_Mag^3*m);
 A(4,12) = (Tmax*lambda_vx*lambda_vy*u)/(lambda_v_Mag^3*m);
 A(4,13) = (Tmax*lambda_vx*lambda_vz*u)/(lambda_v_Mag^3*m);
 A(5,7)  = (Tmax*lambda_vy*u)/(lambda_v_Mag*m^2);
 A(5,11) = (Tmax*lambda_vx*lambda_vy*u)/(lambda_v_Mag^3*m);
 A(5,12) = -(Tmax*u*(lambda_vx^2 + lambda_vz^2))/(lambda_v_Mag^3*m);
 A(5,13) = (Tmax*lambda_vy*lambda_vz*u)/(lambda_v_Mag^3*m);

   A(6,7) = (Tmax*lambda_vz*u)/(lambda_v_Mag*m^2);
  A(6,11) = (Tmax*lambda_vx*lambda_vz*u)/(lambda_v_Mag^3*m);
  A(6,12) = (Tmax*lambda_vy*lambda_vz*u)/(lambda_v_Mag^3*m);
  A(6,13) = -(Tmax*u*(lambda_vx^2 + lambda_vy^2))/(lambda_v_Mag^3*m);

%Factor gravity derivatives; symmetry supplies the other three entries:
       x1 = x+muStar;
       x2 = x1-1;
       b1 = lambda_vx*x1+lambda_vy*y+lambda_vz*z;
       b2 = lambda_vx*x2+lambda_vy*y+lambda_vz*z;
       w1 = 3*(1-muStar)/d5;
       w2 = 3*muStar/r5;
       k1 = 15*(1-muStar)*b1/d7;
       k2 = 15*muStar*b2/r7;
   A(8,1) = -w1*(b1+2*x1*lambda_vx)-w2*(b2+2*x2*lambda_vx)+k1*x1^2+k2*x2^2;
   A(8,2) = -w1*(x1*lambda_vy+y*lambda_vx)-w2*(x2*lambda_vy+y*lambda_vx)+(k1*x1+k2*x2)*y;
   A(8,3) = -w1*(x1*lambda_vz+z*lambda_vx)-w2*(x2*lambda_vz+z*lambda_vx)+(k1*x1+k2*x2)*z;
   A(9,1) = A(8,2);
   A(9,2) = -w1*b1-w2*b2-2*(w1+w2)*y*lambda_vy+(k1+k2)*y^2;
   A(9,3) = -(w1+w2)*(y*lambda_vz+z*lambda_vy)+(k1+k2)*y*z;
  A(10,1) = A(8,3);
  A(10,2) = A(9,3);
  A(10,3) = -w1*b1-w2*b2-2*(w1+w2)*z*lambda_vz+(k1+k2)*z^2;

   A(11,8)  = -1;
   A(12,9)  = -1;
   A(13,10) = -1;
    A(14,7) = (2*Tmax*lambda_v_Mag*u)/m^3;
A(14,11:13) = -Tmax.*lambda_v(:)'.*u./(lambda_v_Mag.*m^2);

 phiDot = A*PHI;
   yDot = [yDot; phiDot(:)];
end

end


