function rv = toPCI(tau,rvND,muStar,tStar,lStar,P)
%% Purpose:
%
%  This routine will take CR3BP dimensionless states and convert them
%  from the rotating barycentric frame of reference to a Planet Centered
%  Inertial Frame of Reference with respect to either primary body
%  as specified by P.  Note that this is not the same as the inertial J2000 
%  system; the planet position is still assumed to be in a circular 
%  restricted orbit with tau = 0 corresponding to the +x axis.
%
%  See 2017. Zimovan. Dissertation for conversions.
%
%% Inputs:
%
%  tau                 [N x 1]              Dimensionless Time
%
%  rvND                [N x 6]              Dimensionless States from
%                                           the CR3BP corresponding to
%                                           tau
%
%  mu                   double              Mass ratio parameter
%                                           mu = m2/(m1+m2)
%
%  tStar                double              Characteristic Time [s]
%
%  lStar                double              Characterisitc Length [km]
%
%  P                    integer             Reference Primary Body 
%                                           1 = Earth Centered
%                                           2 = Moon Centered
%% Outputs:
%
% rv                    [N x 6]             Dimensionalized Position
%                                           and Velocity in a Planet 
%                                           Centered Inertial
%                                           Frame of Reference
%                                           [km,km,km,km/s,km/s,km/s]
%
%% Revision History:
%  Darin C. Koblick                                         (c) 10/21/2025
%  Copyright 2025 Coorbital, Inc.
%% -------------------------- Begin Code Sequence -------------------------
if nargin == 0
  muStar = 0.012150585609624;        %Mass ratio
     tau = linspace(0,2*pi,360)';
   lStar = 389703.264829278;
   tStar = 382981.289129055;
    rvND = repmat([1-muStar,0,0,0,0,0],[size(tau,1) 1]);
      rv = pumpkyn.cr3bp.toPCI(tau,rvND,muStar,tStar,lStar,1);

     figure('Color',[1 1 1]);
     plot3(rv(:,1),rv(:,2),rv(:,3));
     axis equal;
     grid on;
    return;
end

%% Compute the position and velocity of P2:
[rP2,vP2] = pumpkyn.cr3bp.primary2PosVel(tau,muStar);

%% Transform from barycentric to P1:
     rP2 = lStar.*rP2./(1-muStar);
     vP2 = (lStar./tStar).*vP2./(1-muStar);

%% Assemble rotation matrix (3.77 - 3.80)
     xHat = rP2./lStar;                     
     zHat = cross(rP2,vP2,2)./pumpkyn.util.vmag(cross(rP2,vP2,2),2);
     yHat = cross(zHat,xHat,2);

%%  rotation matrix defined using instantaneous rotating axis:
     CR2I = cat(2,permute(xHat,[2 3 1]), ...
                  permute(yHat,[2 3 1]), ...
                  permute(zHat,[2 3 1]));

%% Construct instantaneous angular velocity (3.82):
 thetaDot = pumpkyn.util.vmag(cross(rP2,vP2,2),2)./(lStar.^2);      %rad/s     

%% Translate rotating state such that it's basepoint is located at P2:
rvND(:,1) = rvND(:,1) - (1-muStar);

%% Dimensionalize States:     
         rvND(:,1:3) = rvND(:,1:3).*lStar;                           %km
         rvND(:,4:6) = rvND(:,4:6).*lStar./tStar;                    %km/s

%%  Assemble total transformation matrix:
              CR2I66 = zeros(6,6,size(CR2I,3));
   CR2I66(1:3,1:3,:) =  CR2I;
   CR2I66(4:6,4:6,:) =  CR2I;
   CR2I66(4:6, 1 ,:) = +permute(thetaDot,[3 2 1]).*CR2I(1:3,2,:);
   CR2I66(4:6, 2 ,:) = -permute(thetaDot,[3 2 1]).*CR2I(1:3,1,:);

%% Apply Frame Transformation:
                  rv = pumpkyn.util.multiplyDCM(CR2I66,rvND,2); 

%% Offset from P2 to P1:
if P == 1
           rv(:,1:3) = rv(:,1:3) + rP2;
           rv(:,4:6) = rv(:,4:6) + vP2;
end

end