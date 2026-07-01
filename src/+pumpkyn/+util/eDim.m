function eND = eDim(fND,fSeq)
%% Purpose:
%
%  This routine will reconstruct an N-dimensional matrix from a flattened
%  2-D matrix created by fDim().
%
%% Inputs:
%
%  fND                  [M x K]             Flattened matrix
%
%  fSeq                 struct              Flattening sequence structure
%
%  fSeq.dim             integer             Preserved dimension
%
%  fSeq.postShift       [1 x D]             Matrix size after dimension
%                                           permutation and before reshape
%
%% Outputs:
%
%  eND                  N-D matrix          Reconstructed matrix
%
%% Revision History:
%   Darin Koblick         Created                            (C) 07/19/2012
%   Darin Koblick         Modified                               06/24/2026  
%
%% -------------------------- Begin Code Sequence -------------------------

if isempty(fSeq)
    eND = fND;
    return;
end

dim       = fSeq.dim;
postShift = fSeq.postShift;
     nDim = numel(postShift);
permOrder = localPermOrder(dim,nDim);
      iND = reshape(fND,postShift);
      eND = ipermute(iND,permOrder);

end

function permOrder = localPermOrder(dim,nDim)
    if dim == 0
        permOrder = 1:nDim;
    else
        permOrder = [dim+1:nDim,1:dim];
    end
end