function [x0,tStar,lStar] = fromJ2K(jd0,rv,mu,M,P)
%% Purpose:
%
%  This routine will take J2000 Inertial states with respect to
%  either the Earth or the Moon and convert them to dimensionless
%  states in the rotating barycentric frame such that they can be used
%  with the CR3B equations of motion.
%
%  See 2017. Zimovan. Dissertation.
%
%% Inputs:
%
%  jd0                  [N x 1]             Julian Date Vector (Days)
%
%  rv                   [N x 6]             Dimensionalized Position
%                                           and Velocity in J2000 Inertial
%                                           Frame of Reference
%                                           [km,km,km,km/s,km/s,km/s]
%
%  mu                   double              Mass ratio parameter
%                                           mu = m2/(m1+m2)
%
%  M                    double              Characteristic mass of
%                                           primaries M = m1 + m2
%
%  P                    integer             Reference Primary Body 
%                                           1 = Earth Centered
%                                           2 = Moon Centered
%% Outputs:
%
%  x0                   [N x 6]             Dimensionless States from
%                                           the CR3BP
%
%  tStar                [N x 1]             Characteristic Time (s)
%
%  lStar                [N x 1]             Characteristic Length (km)
%
%% Revision History:
%  Darin C. Koblick                                         (c) 10/01/2025
%  Copyright 2025 Coorbital, Inc.
%% -------------------------- Begin Code Sequence -------------------------
if nargin == 0
      P = 1;
     mu = 0.012150585609624;        %Mass ratio
      M = 5.9736E24 + 7.35E22;      %Characteristic mass (sum of primaries)
    jd0 = juliandate(2020,1,1) + (0:0.1:27.8)';
     x0 = repmat([1-mu,0,0,0,0,0],[size(jd0,1) 1]);
     rv = pumpkyn.cr3bp.toJ2K(jd0,x0,mu,M,P);
     [x0,tStar,lStar] = pumpkyn.cr3bp.fromJ2K(jd0,rv,mu,M,P);
    return;
end
%% Gravitational Constant (km^3/kg/s^2)
  G = 6.67384e-20;             

%% Mass of Both Primaries (kg):
       m2 = mu.*M;                    % mass of P2 (kg)
       m1 = (1-mu).*M;                % mass of P1 (kg)
      mu1 = G.*m1;
      mu2 = G.*m2;

%% Position and Velocity of Moon wrt Earth (405 default model) (km,km/s):
[r12,v12] = planetEphemeris(jd0,'Earth','Moon');

%% Instantaneous Character4istic Quantities:
    lStar = pumpkyn.util.vmag(r12,2);    
    tStar = sqrt((lStar.^3)./(mu1 + mu2));

%% Assemble rotation matrix (3.77 - 3.80)
     xHat = r12./lStar;                     
     zHat = cross(r12,v12,2)./pumpkyn.util.vmag(cross(r12,v12,2),2);
     yHat = cross(zHat,xHat,2);

%%  rotation matrix defined using instantaneous rotating axis:
     CR2I = cat(2,permute(xHat,[2 3 1]), ...
                  permute(yHat,[2 3 1]), ...
                  permute(zHat,[2 3 1]));

%% Construct instantaneous angular velocity (3.82):
 thetaDot = pumpkyn.util.vmag(cross(r12,v12,2),2)./(lStar.^2);     %rad/s     

%%  Assemble total transformation matrix:
              CR2I66 = zeros(6,6,size(CR2I,3));
   CR2I66(1:3,1:3,:) =  CR2I;
   CR2I66(4:6,4:6,:) =  CR2I;
   CR2I66(4:6, 1 ,:) = +permute(thetaDot,[3 2 1]).*CR2I(1:3,2,:);
   CR2I66(4:6, 2 ,:) = -permute(thetaDot,[3 2 1]).*CR2I(1:3,1,:);

%% Transpose the Frame (J2K -> CR3B)
CI2R66 = NaN(size(CR2I66));
for tt=1:size(CR2I66,3)
              CI2R66(:,:,tt) = inv(CR2I66(:,:,tt));
end

%% Offset states to occur at P2:
if P == 1
           rv(:,1:3) = rv(:,1:3) - r12;
           rv(:,4:6) = rv(:,4:6) - v12;
end

%% Apply Frame Transformation:
                  x0 = pumpkyn.util.multiplyDCM(CI2R66,rv,2); 

%% Non-Dimensionalize States:
           x0(:,1:3) = x0(:,1:3)./lStar;
           x0(:,4:6) = x0(:,4:6).*tStar./lStar;

%% Offset states from P2 to rotating barycenter:
              x0(:,1) = x0(:,1) + (1-mu);
end