function phi = elAng(rObs,rTgt,dr,dim3,minElAng)
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
%  minElAng                 Optional scalar elevation threshold (rad), in
%                           [-pi/2, pi/2]. If supplied, return a logical
%                           visibility mask instead of elevation angles.
%
%% Outputs:
%
%  phi                      [N x 1 x M]             Elevation Angle (rad)
%                                                   of target relative
%                                                   to the observer
%                           With minElAng: true at or above the threshold;
%                           false for undefined/nonfinite geometry.
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
        up = rObs-dr;
projection = pumpkyn.util.bsxDot(rObs2Tgt,up,dim3);

if nargin == 5
    if minElAng == 0
        % Only the sign matters at the horizon. Exclude undefined geometry
        % without computing ranges or inverse trigonometric functions.
        valid = all(isfinite(rObs2Tgt),dim3) & all(isfinite(up),dim3) & ...
                any(rObs2Tgt ~= 0,dim3) & any(up ~= 0,dim3);
        phi = valid & projection >= 0;
        return;
    end
end

sinEl = projection ./ (pumpkyn.util.vmag(rObs2Tgt,dim3) .* ...
                       pumpkyn.util.vmag(up,dim3));
% Clamp roundoff at zenith/nadir while preserving undefined angles as NaN.
sinEl(sinEl > 1) = 1;
sinEl(sinEl < -1) = -1;
if nargin == 5
    phi = sinEl >= sin(minElAng);
else
    phi = asin(sinEl);
end
end
