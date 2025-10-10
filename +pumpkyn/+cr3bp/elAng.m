function phi = elAng(rObs,rTgt,dr,dim3)
%% Purpose:
%
%  This routine will determine the elevation angle of rTgt relative to rObs.
%
%% Inputs:
%
%  rObs                       [N x 3]               N number of observers
%
%  rTgt                       [1 x 3 x M]           M number of targets
%                                                   
%
%  dr                       [1 x 3]                 Location of the primary
%                                                   body with respect to
%                                                   the barycenter
%
%  dim3                     integer                 Singleton dimension
%                                                   specifier
%
%% Outputs:
%
%  phi                      [N x 1 x M]             Elevation Angle (rad)
%                                                   of target relative
%                                                   to the observer
%
%% Revision History:
%  Darin C. Koblick                                             08-27-2025
%  Copyright 2025 Coorbital, Inc.
%% ---------------------- Begin Code Sequence -----------------------------
if nargin == 0
                    N = 1500;
                   Np = 15;
                 tau0 = 10*2*pi/9;
 [~,rvTgt,mu,~,lStar] = pumpkyn.cr3bp.getTulip(tau0,Np,-1);
       rP = 1738.1./lStar;
       dr = [1-mu,0,0];
     rObs = pumpkyn.cr3bp.pointSphere(N,dr,rP);
      phi = pumpkyn.cr3bp.elAng(rObs,rvTgt(:,1:3),dr,2);
      idx = phi > 60*pi/180;
      figure('color',[1 1 1]);
      plot3(rObs(:,1),rObs(:,2),rObs(:,3),'.k'); hold on;
      plot3(rvTgt(:,1),rvTgt(:,2),rvTgt(:,3),'.b');
      plot3(rObs(idx,1),rObs(idx,2),rObs(idx,3),'.r');
      axis equal;
    return;
end
rObs2Tgt = rTgt-rObs;
   theta = pumpkyn.util.bsxAng(rObs2Tgt,rObs-dr,dim3);   %Zenith Angle
     phi = (90 - theta).*pi/180;            %90 deg - Zenith = Elevation
end