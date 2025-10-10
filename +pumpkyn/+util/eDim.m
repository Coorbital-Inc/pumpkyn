function eND = eDim(fND,fSeq)
%% Purpose:
%Take any flattened 2-D matrix in MATLAB and convert it back into its
%original multi-dimensional form before flatening. This is often times 
%necessary when writing complex operations on multi-dimensional matrices.  
%It is also desired that after flattening, the dimension that is preserved 
%has the correct sequence. This is especially important for vector 
%processing. Once flattened, and an operation has been performed on the 2-D
%matrix, often times the 2-D matrix will need to be converted back to the
%original multi-dimensional matrix.  This function will extract the
%multi-dimensional matrix from the flattened 2-D matrix using fDim().
%
%% Inputs:
% fND           [N x size(ND,dim)]                  A 2-D matrix with
%                                                   N corresponding to the
%                                                   same number of elements
%                                                   in the ND matrix but
%                                                   with a single dimension
%                                                   preserved.
%
% fSeq           struct                             fSeq.dim is the initial
%                                                   location of the desired
%                                                   dimension(s) before the
%                                                   shift sequence.
%                                                   
%                                                   fSeq.postShift is the
%                                                   matrix dimensions
%                                                   before the reshaping
%                                                   sequence.
%
%
%% Outputs:
% eND            [O x P x Q x R x S x ...]          Unflattened matrix with 
%                                                   any number and order of
%                                                   dimensions.
%
%
%% Created By Darin Koblick (C) 07/19/2012
iND = reshape(fND,fSeq.postShift);
%inject a singleton dimension:
if numel(size(iND)) ~= numel(fSeq.postShift)
    %shift dimensions to the right and add singletons along the way:
    eND = shiftdim(iND,numel(size(iND))-numel(fSeq.postShift));
    %Add a circular shift to the left and wrap the leading edges:
    eND = shiftdim(eND,numel(fSeq.postShift)-(fSeq.dim)+1);
else
    %Shift dimensions to the left and wrap the leading edges:
    eND = shiftdim(iND,numel(fSeq.postShift)-(fSeq.dim));
end
end