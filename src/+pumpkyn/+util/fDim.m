function [fND,fSeq] = fDim(ND,dim)
%% Purpose:
%
%  This routine will flatten an N-dimensional matrix into a 2-D matrix
%  while preserving the specified dimension as the column dimension.
%
%% Inputs:
%
%  ND                   N-D matrix          Input matrix
%
%  dim                  integer             Dimension to preserve as columns
%                                           dim = 0 flattens to [numel x 1]
%
%% Outputs:
%
%  fND                  [M x size(ND,dim)]  Flattened matrix
%
%  fSeq                 struct              Flattening sequence structure
%
%  fSeq.dim             integer             Preserved dimension
%
%  fSeq.postShift       [1 x D]             Matrix size after dimension
%                                           permutation and before reshape
%
%% Revision History:
%   Darin Koblick         Created                            (C) 07/19/2012
%   Darin Koblick         Modified                               06/24/2026  
%
%% -------------------------- Begin Code Sequence -------------------------

if isempty(ND)
    fND  = ND;
    fSeq = [];
    return;
end

if nargin < 2 || isempty(dim)
    dim = 0;
end

if ~isscalar(dim) || dim < 0 || dim ~= floor(dim)
    error('fDim:BadDim','dim must be a nonnegative integer scalar.');
end

% Preserve original behavior: invalid dim means flatten to column
if dim > ndims(ND)
    dim = 0;
end

nDim = ndims(ND);

sizeND = size(ND);
sizeND(end+1:nDim) = 1;

permOrder = localPermOrder(dim,nDim);

ND = reshape(ND,sizeND);
iND = permute(ND,permOrder);

fSeq.dim       = dim;
fSeq.postShift = sizeND(permOrder);

if dim == 0
    fND = reshape(iND,numel(ND),1);
else
    fND = reshape(iND,numel(ND)/sizeND(dim),sizeND(dim));
end

end

function permOrder = localPermOrder(dim,nDim)

if dim == 0
    permOrder = 1:nDim;
else
    permOrder = [dim+1:nDim,1:dim];
end

end