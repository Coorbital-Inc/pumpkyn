function [s,cs] = arclength(r,dim3)
%% Purpose:
%
%  This routine computes the total length of a line given in vector form 
%  as [N x M] where N is the number of points, and M is the singleton 
%  dimension 
%
%% Inputs:
%
%  r                    [N x M]                     N-D Vectors
%
%  dim3                 integer                     Singleton dimension
%                                                   describing location of
%                                                   M in the vector array
%
%
%% Outputs:
%
%  s                    double                      Total Arclength of line
%
% cs                    [N x 1]                     Cummulative Arclength of
%                                                   the line
%
%% Revision History:
%  Darin C. Koblick                                         (c) 02-04-2020
%% --------------------------- Begin Code Sequence ------------------------
if nargin == 0
     r = [1 0 0
          2 0 0
          3 0 0
          4 0 0];
      
     [s,cs] = arclength(r,2); 
     return;
end
       r = pumpkyn.util.fDim(r,dim3);
  deltaR = diff(r,1,1);
      cs = [0; cumsum(sqrt(sum(deltaR.^2,2)),1)];
       s = cs(end);
end