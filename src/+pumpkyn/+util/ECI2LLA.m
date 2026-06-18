function LLA = ECI2LLA(JD,r_ECI,dim)
% Purpose:                                                                
% Convert ECI (CIS, Epoch J2000.0) Coordinates to Latitude/Longitude and  
% altitude                                                                
%
% Inputs:                                                                 
%-------                                                                  
%JD                     [N x 1]                         Julian Date Vector
%
%r_ECI                  [N x 3]                         Position Vector
%                                                       in ECI coordinate
%                                                       frame of reference
%
%dim                    int                             singleton dimension 
%                                                       specifier
% Outputs:
%---------                                                                %
% LLA                   [N x 3]                         Latitude,Longitude,
%                                                       and Altitude
%                                                       in same dimensions
%                                                       as the ECI position
%                                                       vectors. (rad,rad,
%                                                       and km
%                                                       respectively)
%
% References:
%-------------
%Orbital Mechanics with Numerit, http://www.cdeagle.com/omnum/pdf/csystems.pdf
%
%
% Function Dependencies:
%------------------
% JD2GMST
%------------------------------------------------------------------       %
% Programed by Darin Koblick  01-08-2014                                  %
%----------------------- Begin Code Sequence -----------------------------%

%Flatten the input matricies
mSize = size(r_ECI);
mSize(dim) = 1;
JD = bsxfun(@times,JD,ones(mSize));
JD = pumpkyn.util.fDim(JD,dim);
[r_ECI,fSeq] = pumpkyn.util.fDim(r_ECI,dim);
%Call the routine to convert ECI to ECEF:
r_ECEF = pumpkyn.util.ECItoECEF(JD,r_ECI,r_ECI,r_ECI,2);
%Find the lat/lon/alt:
LLA = pumpkyn.util.ECEF2LLA(r_ECEF,2);
%Convert back to original matrix dimensions:
LLA = pumpkyn.util.eDim(LLA,fSeq);
end