function [r,v] = LLA2ECI(jd,lla,dim3)
%% Purpose:
%LLA2ECI - convert latitude, longitude, and altitude to
%            earth-centered, earth-inertial (ECI) cartesian
% 
% USAGE:
% [r,v] = LLA2ECI(jd, [lat,lon,alt],dim3)
%
%
%% Inputs:
% jd                [M x N x O x ... 1]         Julian Date Vector
%                                               correpsonding to each
%                                               element contained in the
%                                               latitude, longitude and
%                                               altitude matrix
%
% lla               [M x N x O x ... 3]          Multi dimensional array
%                                                where the dim3
%                                                specification denotes the
%                                                following terms:
%                                                geodetic latitude (radians), 
%                                                longitude (radians), 
%                                                height above spherical 
%                                                earth (km)
%
% dim3              double                      Dimension where lat,lon,alt
%                                               can be found
% 
%
%% Outputs:
% r              [M x N x O x ... 3]            Multi dimensional position
%                                               array where dim3 
%                                               specification denotes the 
%                                               following terms:
%                                               ECI X-coordinate (km)
%                                               ECI Y-coordinate (km)
%                                               ECI Z-coordinate (km)
%
% v             [M x N x O x ... 3]             Multi dimensional velocity
%                                               array where dim3 
%                                               specification denotes the 
%                                               following terms:
%                                               ECI Xdot-coordinate (km/s)
%                                               ECI Ydot-coordinate (km/s)
%                                               ECI Zdot-coordinate (km/s)
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
%%  Revision History:
%  Darin C. Koblick                                         (c) 10-24-2014
%% ------------------- Begin Code Sequence --------------------------------

%% Flatten dimensions:
 [lla,fSeq] = pumpkyn.util.fDim(lla,dim3);
         jd = pumpkyn.util.fDim(jd,dim3);
%% Step One: Convert LLA to ECEF:
r_ECEF = pumpkyn.util.LLA2ECEF(lla,2)';

%% Step Two: Convert ECF to ECI:
[r,v] = pumpkyn.util.ECEFtoECI(jd(:)',r_ECEF,zeros(size(r_ECEF)),zeros(size(r_ECEF)));

%% Convert the position and velocity back to original dimension specs:
r = pumpkyn.util.eDim(r',fSeq);
v = pumpkyn.util.eDim(v',fSeq);

end