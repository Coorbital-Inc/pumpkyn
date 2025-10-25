function out = vmag(in,dim)
%% Purpose:
% Quickly compute the vector magnitude of any given N-D Vector. This is
% equilivent of finding the 2-norm function in MATLAB.
%
%% Inputs:
%  in               [N x M x O x ... P]                      Any dimensions
%
%  dim              int                                      Dimension
%                                                            specification
%                                                            of X,Y,Z
%
%% Outputs:
% out              [N x M x O x .... 1]                     Resulting 
%                                                           Vector
%                                                           magnitude 
%% Revision History:
% Darin C. Koblick                                          (c) 07-17-2013 
%% --------------------- Begin Code Sequence ------------------------------
out = sqrt(sum(in.^2,dim));
end