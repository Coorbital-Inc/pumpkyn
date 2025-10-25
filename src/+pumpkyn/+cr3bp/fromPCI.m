function rvND = fromPCI(t,rv,muStar,tStar,lStar,P)
%% Purpose:
%
%  This routine will take Planet Ceneterd Inertial states with respect to
%  either Primary and convert them to dimensionless
%  states in the rotating barycentric frame such that they can be used
%  with the CR3B equations of motion. This is a complimentary function
%  to the toPCI.m routine, not to be used or confused with fromJ2K which
%  is epoch dependent and transforms to the Moon Centered or Earth Centered
%  inertial frame.
%
%  See 2017. Zimovan. Dissertation for conversions.
%
%% Inputs:
%
%   t                   [N x 1]             Dimensionalized Time Vector [s]
%
%  rv                   [N x 6]             Dimensionalized Position
%                                           and Velocity in PCI                                           Frame of Reference
%                                           [km,km,km,km/s,km/s,km/s]
%
%  muStar               double              Mass ratio parameter
%                                           mu = m2/(m1+m2)
%
%  tStar                double              Characteristic Time [s]
%
%  lStar                double              Characterisitc Length [km]
%
%
%  P                    integer             Reference Primary Body 
%                                           1 = Earth Centered
%                                           2 = Moon Centered
%
%% Outputs:
%
%  x0                   [N x 6]             Dimensionless States from
%                                           the CR3BP
%
%% Revision History:
%  Darin C. Koblick                                         (c) 10/01/2025
%  Copyright 2025 Coorbital, Inc.
%% -------------------------- Begin Code Sequence -------------------------
if nargin == 0
  muStar = 0.012150585609624;        %Mass ratio
     tau = linspace(0,2*pi,360)';
   lStar = 389703.264829278;
   tStar = 382981.289129055;
    rvND = repmat([1-muStar,0,0,0,0,0],[size(tau,1) 1]);
       P = 1;
      rv = pumpkyn.cr3bp.toPCI(tau,rvND,muStar,tStar,lStar,P);
  rvND_2 = pumpkyn.cr3bp.fromPCI(tau.*tStar,rv,muStar,tStar,lStar,P);
  figure('Color',[1 1 1]);
  plot(pumpkyn.util.vmag(rvND_2-rvND,2));
  grid on;
  ylabel('RSS Error');
    return;
end

%% Compute the position and velocity of P2:
[rP2,vP2] = pumpkyn.cr3bp.primary2PosVel(t./tStar,muStar);

%% Transform from barycentric to P1:
      rP2 = lStar.*rP2./(1-muStar);
      vP2 = (lStar./tStar).*vP2./(1-muStar);

%% Construct instantaneous angular velocity (3.82):
 thetaDot = pumpkyn.util.vmag(cross(rP2,vP2,2),2)./(lStar.^2);      %rad/s

%% Assemble rotation matrix (3.77 - 3.80)
     xHat = rP2./lStar;                     
     zHat = cross(rP2,vP2,2)./pumpkyn.util.vmag(cross(rP2,vP2,2),2);
     yHat = cross(zHat,xHat,2);

%%  rotation matrix defined using instantaneous rotating axis:
     CR2I = cat(2,permute(xHat,[2 3 1]), ...
                  permute(yHat,[2 3 1]), ...
                  permute(zHat,[2 3 1]));
     
%%  Assemble total transformation matrix:
              CR2I66 = zeros(6,6,size(CR2I,3));
   CR2I66(1:3,1:3,:) =  CR2I;
   CR2I66(4:6,4:6,:) =  CR2I;
   CR2I66(4:6, 1 ,:) = +permute(thetaDot,[3 2 1]).*CR2I(1:3,2,:);
   CR2I66(4:6, 2 ,:) = -permute(thetaDot,[3 2 1]).*CR2I(1:3,1,:);

%% Transpose the Frame (PCI -> CR3B)
CI2R66 = NaN(size(CR2I66));
for tt=1:size(CR2I66,3)
              CI2R66(:,:,tt) = inv(CR2I66(:,:,tt));
end

%% Offset states to occur at the P2:
if P == 1
           rv(:,1:3) = rv(:,1:3) - rP2;
           rv(:,4:6) = rv(:,4:6) - vP2;
end

%% Apply Frame Transformation:
                 rv = pumpkyn.util.multiplyDCM(CI2R66,rv,2); 

%% Non-Dimensionalize States:
         rvND(:,1:3) = rv(:,1:3)./lStar;
         rvND(:,4:6) = rv(:,4:6).*(tStar./lStar);

%% Offset states from P2 to rotating barycenter:
           rvND(:,1) = rvND(:,1) + (1-muStar);
end