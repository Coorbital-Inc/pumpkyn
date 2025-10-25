function [rvDot,aP1,aP2] = threeBody(t,rv,muStar,M,tStar,lStar,P,dim3)
%% Purpose:
%
%  This routine will propagate a satellite using the cartesian representation
%  of the two-body equations of motion with an additional perturabtion
%  term representing the third body.  All inputs and outputs are 
%  dimensionalized and convertable to/from CR3BP via fromPCI() and toPCI().
%  
%  Note that this assumes an epoch such that there is an x-axis crossing
%  when t=0, and that the second primary, P2, is in a circular non-inclined
%  orbit (Circular Restricted).
%
%% Inputs:
%
%  t                        [N x 1]             Dimensionalized Time
%
%  rv                       [N x 6]             Dimensionalized position
%                                               and velocity (km, and km/s)
%
%  muStar                    double             Gravitaional Mass Ratio of
%                                               two primaries such that
%                                               muStar = m2/(m1+m2)
%
%  M                        double              Characteristic mass of
%                                               primaries M = m1 + m2
%
%  tStar                    double              Characteristic Time
%
%
%  lStar                    double              Characteristic Length
%
%
%  P                        integer             Primary Body of coorinate
%                                               origin:
%                                               1 = P1 Centered
%                                               2 = P2 Centered
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
%  Darin C. Koblick                                         (c) 10-25-2025
%  Copyright 2025 Coorbital, Inc.
%% -------------------------- Begin Code Sequence -------------------------
if nargin == 0
        M = 5.9736E24 + 7.35E22;
        P = 2;
       Np = 7;     %Number of petals
     tau0 = 2*pi; 
       pm = +1;
   [tau0,rvND0,muStar,tStar,lStar] = pumpkyn.cr3bp.getTulip(tau0,Np,pm);
      rv0 = pumpkyn.cr3bp.toPCI(0,rvND0,muStar,tStar,lStar,P);
     PHI0 = eye(6);
     opts = odeset('absTol',1e-13,'relTol',1e-13);
   [t,rv] = ode113(@pumpkyn.cr3bp.threeBody,[0 tau0*tStar],[rv0,PHI0(:)'],opts,muStar,M,tStar,lStar,P,1);

   %Convert back to the CR3BP frame:
     rvND = pumpkyn.cr3bp.fromPCI(t,rv(:,1:6),muStar,tStar,lStar,P);

   figure('Color',[1 1 1]);
   plot3(rvND(:,1),rvND(:,2),rvND(:,3));
   grid on;
   axis equal;

   return;
end

%% Flatten inputs:
   [rv,fSeq] =  pumpkyn.util.fDim(rv,dim3);
       rvDot =  NaN(size(rv));

%% Gravitational Constant (km^3/kg/s^2)
        G = 6.67384e-20;             
      tau = t./tStar;

%% Mass of Both Primaries (kg):
       m2 = muStar.*M;                    % mass of P2 (kg)
       m1 = (1-muStar).*M;                % mass of P1 (kg)
      mu1 = G.*m1;
      mu2 = G.*m2;

%% Determine P2 Position wrt P1:
      rP2 = pumpkyn.cr3bp.primary2PosVel(tau,muStar);
      rP2 = lStar.*rP2./(1-muStar);                   %km

%% Ensure that we are P1 or P2 centered:
if P == 2
%P2 Centered, use P2 as (0,0,0)
          rP2 = -rP2;
    [mu1,mu2] = deal(mu2,mu1);
end

%% Three-Body Equations of Motion (Vallado 8.6.3 Eq. 8-34):
rvDot(:,1:3) =  rv(:,4:6);
      rP1Sat =  rv(:,1:3);
      rSatP2 = rP2 - rP1Sat; 
   rP1SatMag = sqrt(dot(rP1Sat,rP1Sat,2));
   rSatP2Mag = sqrt(dot(rSatP2,rSatP2,2));
      rP2Mag = sqrt(dot(rP2,rP2,2));
         aP1 = -mu1.*rP1Sat./rP1SatMag.^3;
         aP2 =  mu2.*(rSatP2./rSatP2Mag.^3 - rP2./rP2Mag.^3);
rvDot(:,4:6) = aP1 + aP1;

%% Assemble the STM Derivative:
if size(rv,2) >= 42
    %Extract the State Transition Matrix:
      PHI = reshape(rv(:,7:42)',[6 6]);
     % Jacobian of vector field:
     % d_rvDot/d_rv
     dadr = compute_dadr(rP1Sat,rP2,mu1,mu2);
        A = [zeros(3),eye(3); dadr, zeros(3)];
   PHIdot = A*PHI;
    rvDot = [rvDot(:,1:6), PHIdot(:)'];
end

%% Package up outputs:
       rvDot = pumpkyn.util.eDim(rvDot,fSeq);
end

%d_a/d_r
function dadr = compute_dadr(r,rP2,mu1,mu2)
% syms x y z xdot ydot zdot xP2 yP2 zP2 mu1 mu2 real;
% 
%          rP2 = [xP2 yP2 zP2];
%            r =  [x, y, z];
%       rSatP2 = rP2 - r; 
%    rP1SatMag = sqrt(dot(r,r,2));
%    rSatP2Mag = sqrt(dot(rSatP2,rSatP2,2));
%       rP2Mag = sqrt(dot(rP2,rP2,2));
%            a = -mu1.*r./rP1SatMag.^3 + ...
%                 mu2.*(rSatP2./rSatP2Mag.^3 - rP2./rP2Mag.^3);
%            v = [xdot,ydot,zdot];
% 
%  %Jacobian of vector field:
%  % d_va/d_rv
%  f = [v,a];
% 
%  % Compute the Jacobian matrix of the vector field
%     J = jacobian(f, [r,v]);
%  dadr = simplify(J(4:6,1:3));

        x = r(:,1);
        y = r(:,2);
        z = r(:,3);
      xP2 = rP2(:,1);
      yP2 = rP2(:,2);
      zP2 = rP2(:,3);
     dadr = NaN(3,3);

dadr(1,1) = (3*mu1*x^2)/(x^2 + y^2 + z^2)^(5/2) - ...
            (mu2*(- 2*x^2 + 4*x*xP2 - 2*xP2^2 + y^2 - ...
            2*y*yP2 + yP2^2 + z^2 - 2*z*zP2 + zP2^2)) ...
            /((x - xP2)^2 + (y - yP2)^2 + (z - zP2)^2)^(5/2) -  ...
            mu1/(x^2 + y^2 + z^2)^(3/2);

dadr(1,2) = (3*mu1*x*y)/(x^2 + y^2 + z^2)^(5/2) + ...
            (3*mu2*(2*y - 2*yP2)*(x - xP2))/(2*((x - xP2)^2 + ...
            (y - yP2)^2 + (z - zP2)^2)^(5/2));

dadr(1,3) = (3*mu1*x*z)/(x^2 + y^2 + z^2)^(5/2) + ...
            (3*mu2*(2*z - 2*zP2)*(x - xP2))/ ...
            (2*((x - xP2)^2 + (y - yP2)^2 + (z - zP2)^2)^(5/2));

dadr(2,1) = dadr(1,2);

dadr(2,2) = (3*mu1*y^2)/(x^2 + y^2 + z^2)^(5/2) - ...
            (mu2*(x^2 - 2*x*xP2 + xP2^2 - 2*y^2 + 4*y*yP2 - ...
             2*yP2^2 + z^2 - 2*z*zP2 + zP2^2))/ ...
             ((x - xP2)^2 + (y - yP2)^2 + (z - zP2)^2)^(5/2) - ...
             mu1/(x^2 + y^2 + z^2)^(3/2);

dadr(2,3) = (3*mu1*y*z)/(x^2 + y^2 + z^2)^(5/2) + ...
            (3*mu2*(2*z - 2*zP2)*(y - yP2))/ ...
            (2*((x - xP2)^2 + (y - yP2)^2 + (z - zP2)^2)^(5/2));

dadr(3,1) = dadr(1,3);

dadr(3,2) = dadr(2,3);

dadr(3,3) = (3*mu1*z^2)/(x^2 + y^2 + z^2)^(5/2) - ...
            (mu2*(x^2 - 2*x*xP2 + xP2^2 + y^2 - 2*y*yP2 + yP2^2 - 2*z^2 + ...
            4*z*zP2 - 2*zP2^2))/((x - xP2)^2 + (y - yP2)^2 + ...
            (z - zP2)^2)^(5/2) - mu1/(x^2 + y^2 + z^2)^(3/2);
end