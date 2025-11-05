function [sND1,sND2] = sDim(ND,dim)
%% Purpose:
%
%  Given a single input matrix with it's singleton specifier defined by
%  dim, this routine will split it in half with the first half output as
%  sND1, and the second half output as sND2
%
%% Inputs:
%
% ND            [O x P x Q x R x S x ...]           Any matrix with any
%                                                   number, N, of
%                                                   dimensions.
%
% dim           double                              Specify a single
%                                                   dimension in which to
%                                                   preserve when creating
%                                                   the columns of the 2-D
%                                                   flattened matrix.
%
%% Outputs:
%
%  sND1          [O x P x Q x R x S x ...]        First Split matrix with
%                                                 number, N, dimensions
%                                                 along the singleton dim 
%
%  sND2          [O x P x Q x R x S x ...]        Second Split matrix with 
%                                                 N dimensions along the 
%                                                 singleton dim 
%
%% Revision History:
%  Darin C. Koblick                                         (c) 09-21-2020
%% ------------------------ Begin Code Sequence ---------------------------

if nargin == 0
            ND = repmat([1 2 3 4 5 6],[10 1]); 
   [sND1,sND2] = pumpkyn.util.sDim(ND,2);
   return;
end
          [ND,fSeq] = pumpkyn.util.fDim(ND,dim);       % Flattened
               sND1 = ND(:,1:end/2);      % First half
               sND2 = ND(:,end/2+1:end);  % Second half
fSeq.postShift(end) = fSeq.postShift(end)/2;
               sND1 = pumpkyn.util.eDim(sND1,fSeq);
               sND2 = pumpkyn.util.eDim(sND2,fSeq);
end