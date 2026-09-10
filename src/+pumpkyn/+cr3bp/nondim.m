function [tau0,x0] = nondim(tau0,x0,tStar,lStar,dim3)
%% Purpose:
%
%  This routine will take dimensional cr3bp states in the rotating 
%  barycentric frame and non-dimensionalize them according to the user provided
%  characteristic time and length values respectively.
%
%% Inputs:
%
%  tau0                 [N x 1]             Dimensional Period
%
%
%  x0                   N-D array           Dimensional [x,y,z,vx,vy,vz]
%                                           along dim3. With 42 or more
%                                           components, 7:42 contain Phi(:)
%                                           in MATLAB column-major order.
%
%  tStar                positive scalar     Characteristic Time (e.g. s)
%
%  lStar                positive scalar     Characteristic Length (e.g. km)
%
%  dim3                 integer             State-component dimension:
%                                           2 for N x 6 or N x 42;
%                                           1 for a column state vector.
%
%% Outputs:
%
%  tau0                 [N x 1]             Dimensionless Period
%
%
%  x0                   N-D array           Dimensionless position, velocity
%                                           and STM (when present). Other
%                                           components remain unchanged.
%
%% Revision History:
%  Darin C. Koblick                                         (c) 09/10/2026
%  Copyright 2026 Coorbital, Inc.
%% -------------------------- Begin Code Sequence -------------------------
if nargin == 0
    tStar = 382981.289129055;
    lStar = 389703.264829278;
     tau0 = 2743458.2027111659;
       x0 = [446086.10905823804	
            -9.9655787903683696e-15	
             35714.428021262051	
             0.0092656172279980773	
            -0.35454425599683487	
             0.30745467347379113]';
 [tau0,x0] = pumpkyn.cr3bp.nondim(tau0,x0,tStar,lStar,2);
 return;
end
%Non-dimensionalize time while preserving its input shape:
     tau0 = tau0./tStar;
%Flatten the input dimensions for x0:
[x0,fSeq] = pumpkyn.util.fDim(x0,dim3);
%Scale position and velocity separately:
x0(:,1:3) = x0(:,1:3)./lStar;
x0(:,4:6) = x0(:,4:6).*(tStar/lStar);
%Scale the two off-diagonal STM blocks directly in packed Phi(:) order:
if size(x0,2) >= 42
    x0(:,[25:27 31:33 37:39]) = x0(:,[25:27 31:33 37:39])./tStar;
    x0(:,[10:12 16:18 22:24]) = x0(:,[10:12 16:18 22:24]).*tStar;
end
%Reshape the outputs:
       x0 = pumpkyn.util.eDim(x0,fSeq);
end
