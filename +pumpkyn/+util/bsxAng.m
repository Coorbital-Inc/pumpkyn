function ang = bsxAng(A,B,dim3)
%% Purpose:
% Take two multi-dimensional arrays which contain vector data and compute 
% the angle between each vector.
%
%% Inputs:
%  A                     [N x M x O x ... ]               Multi-Dimensional
%                                                         vector (x,y,z)
%                                                         3-D
%
%  B                     [N x M x O x ... ]               Multi-Dimensional
%                                                         vector (x,y,z)
%                                                         3-D
%
%  dim3                   int                             Specificaion of
%                                                         the x,y,z
%                                                         dimension along
%                                                         vectors A and B
%
%% Outputs:
%  ang                  [N x M x O x ... 1]               Angle between
%                                                         vectors A and B
%                                                         (deg)
%% Revision History:
%  Darin C. Koblick                                               (c) 2013
%% -------------------- Begin Code Sequence -------------------------------

%Compute the Angle Between arrays:
ang = real(acosd(bsxfun(@rdivide,pumpkyn.util.bsxDot(A,B,dim3),bsxfun(@times,pumpkyn.util.vmag(A,dim3),pumpkyn.util.vmag(B,dim3)))));

end