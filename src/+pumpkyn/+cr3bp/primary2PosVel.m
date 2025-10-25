function [r,v] = primary2PosVel(tau,muStar)
%% Purpose:
%
%  This routine will compute the dimensionless position and velocity of 
%  the second primary assuming the circular restricted three-body 
%  assumption such that its orbital period is equal to 2*pi.
%
%  tau = 0 corresponds to the x-axis, such that, r = [1,0,0]
%
%% Inputs:
%
%  tau                  [N x 1]             Dimensionless Time from
%                                           tau = 0
%
%  muStar               double              Dimensionless Mass Ratio of
%                                           the two primaries
%                                           muStar = m2/(m2+m1)
%
%% Outputs:
%
%  r                    [N x 3]             Dimensionless Position of P2
%                                           corresponding to tau with 
%                                           respect to barycenter
%
%  v                    [N x 3]             Dimensionless Velocity of P2
%                                           corresponding to tau with
%                                           respect to barycenter
%% Revision History:
%  Darin C. Koblick                                         (c) 10/21/2025
%  Copyright 2025 Coorbital, Inc.
%% -------------------------- Begin Code Sequence -------------------------
if nargin == 0
     tau = linspace(0,pi,360)';
  muStar = 0.012150585609624;        %Mass ratio
   [r,v] = pumpkyn.cr3bp.primary2PosVel(tau,muStar);
   return;
end
%Assume we have + orbit momentum:
r = [+cos(tau), sin(tau), tau.*0].*(1-muStar);
%Temporal Derivative of position vector:
v = [-sin(tau), cos(tau), tau.*0].*(1-muStar);
end