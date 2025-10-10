function maxGapTimes = maxGaps(tau,r,rPts,rPt0,minElAng,Ns,dTau)
%% Purpose:
%
%  Computes the maximum surface gap times for a set of discretized points 
%  representing P2, provided by rPts as seen by the satellites in an 
%  orbital configuration denoted by r.
%
%
%% Inputs:
%
%  tau                      [N x 1]                 Dimensionless Time
%
%  r                        [N x 3 x M]             Dimensionless Position
%                                                   Vector of all satellites
%                                                   with respect to tau
%
%  rPts                     [P x 3]                 Dimensionless Pos
%                                                   of all surface points
%                                                   on primary body
%
%  rPt0                     [1 x 3]                 Dimensionless Pos
%                                                   of Primary Body
%
%  minElAng                 double                  Minimum Acceptable
%                                                   Elevation Angle (rad)
%
%  Ns                       integer                 Minumum Number of
%                                                   satellites in view
%                                                   at one time to be
%                                                   considered an
%                                                   observation.
%
%  dTau                     double                  Dimensionless time
%                                                   step to discretize
%                                                   the satellite state
%                                                   vector (lower values
%                                                   equal higher
%                                                   resolution)
%
%% Outputs:
%
%  maxGapTimes             [P x 1]                  Maximum Time
%                                                   between observations
%                                                   of at least N satellite
%                                                   and a specific surface
%                                                   point (dimensionless)
%
%% Revision History:
%  Darin C. Koblick                                         (c) 10-03-2025
%  Copyright 2025 Coorbital, Inc.
%% ----------------------- Begin Code Sequence ----------------------------

%% Interpolate based on discretized timesteps:
 tau_i = (0:dTau:tau(end))';
     r = interp1(tau,r,tau_i,'spline');
   tau = tau_i;

%% Compute Elevation Angle:
 elAng = pumpkyn.cr3bp.elAng(permute(rPts,[4 2 3 1]),r,rPt0,2);

%% For each surface point, compute the gaps:
    maxGapTimes = zeros(size(rPts,1),1);
for tp=1:size(rPts,1)
         satAngs = permute(elAng(:,1,:,tp),[1 3 2]);
     satsPerStep = sum(satAngs >= minElAng,2);
        obsTimes = [0; tau(satsPerStep >= Ns); tau(end)];
        gapTimes = diff(obsTimes);
       gapTimes(gapTimes < dTau*1.5) = [];        %Remove gaps near dTau
       if ~isempty(gapTimes)
         maxGapTimes(tp,1) = max(gapTimes);
       end
end

end