function [xyz] = LLA2ECEF(lla,dim3)
%% Purpose:
%LLA2ECEF - convert latitude, longitude, and altitude to
%            earth-centered, earth-fixed (ECEF) cartesian
% 
% USAGE:
% [x,y,z] = lla2ecef(lat,lon,alt)
%
%
%% Inputs:
% lla               [M x N x O x ... 3]          Multi dimensional array
%                                                where the dim3
%                                                specification denotes the
%                                                following terms:
%                                                geodetic latitude (radians), 
%                                                longitude (radians), 
%                                                height above Earth's Surface
%                                                (km)
%
% dim3              double                      Dimension where lat,lon,alt
%                                               can be found
% 
%
%% Outputs:
% xyz              [M x N x O x ... 3]          Multi dimensional array
%                                               where dim3 specification
%                                               denotes the following
%                                               terms:
%                                               ECEF X-coordinate (km)
%                                               ECEF Y-coordinate (km)
%                                               ECEF Z-coordinate (km)
%
%
% Notes: This function assumes the WGS84 model.
%        Latitude is customary geodetic (not geocentric).
% 
% Source: "Department of Defense World Geodetic System 1984"
%         Page 4-4
%         National Imagery and Mapping Agency
%         Last updated June, 2004
%         NIMA TR8350.2
% 
% Michael Kleder, July 2005
% Modified by Darin Koblick to support multi-dimensional inputs
%             June 2013
%% Begin Code Sequence
if nargin == 1
   dim3 = 2; 
end
%Flatten the input matrix
[lla,Seq] = pumpkyn.util.fDim(lla,dim3);
% WGS84 ellipsoid constants:
a = 6378.137;
e = 8.1819190842622e-2;
% intermediate calculation
% (prime vertical radius of curvature)
N = a./sqrt(1 - e^2 .* sin(lla(:,1)).^2);
xyz = cat(2,(N+lla(:,3)) .* cos(lla(:,1)) .* cos(lla(:,2)), ...
         (N+lla(:,3)) .* cos(lla(:,1)) .* sin(lla(:,2)), ...
         ((1-e^2) .* N + lla(:,3)) .* sin(lla(:,1)));
%Reconstruct the output matrix into the same shape as the input matrix
xyz = pumpkyn.util.eDim(xyz,Seq);
return