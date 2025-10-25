function C = bsxDot(A,B,dim3)
%% Purpose:
% Take two multi-dimensional arrays which contain vector data and compute 
% the dot product between each vector.
%
%% Inputs:
%  A                     [N x M x O x ... 3]              Multi-Dimensional
%                                                         vector (x,y,z)
%                                                         3-D
%
%  B                     [N x M x O x ... 3]              Multi-Dimensional
%                                                         vector (x,z,z)
%                                                         3-D
%
%  dim3                   int                             Specificaion of
%                                                         the x,y,z
%                                                         dimension along
%                                                         vectors A and B
%
%% Outputs:
%  C                  [N x M x O x ... 1]                 Dot product
%                                                         between vectors A
%                                                         and B
%% Revision History:
%  Darin C. Koblick                                       (c) 06-19-2015
%% -------------------- Begin Code Sequence ------------------------------- 
if nargin == 0
   A = rand([3 10 20 40 50 10]);
   B = rand([3 10 20 40 50 10]);
   tic;
   C = bsxDot(A,B,4);
   disp(num2str(toc));
   tic;
   D = dot(A,B,4);
   disp(num2str(toc));
   isequal(C,D)
   C = [];
   return;
end

if ~exist('dim3','var')
   dim3 = 2; 
end
    C = sum(A.*B,dim3);
end