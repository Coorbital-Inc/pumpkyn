function s = multiplyDCM(DCM,r,dim3)
%% Purpose:
%  Take a L x M x N DCM and multiply it by a N x M vector
%  The output is a vector will be of dimension N x M
%
%% Inputs:
%  DCM                  [L x M x N]            Direction Cosine Matrix
%
%  r                    [N x M]                Initial Vector in which
%                                              you would like to
%                                              multiply by the DCM
%
% dim3                  int                    Singleton dimension
%                                              specifier of the input
%                                              vector, r. (e.g. M)
%
%% Outputs:
%  s                    [N x L]                Vector r rotated by the
%                                              DCM, that is: s = DCM*r
%% Revision History:
%  Darin C. Koblick                                         (c) 01-30-2023
%% --------------------- Begin Code Sequence ------------------------------
if nargin == 0
    DCM = repmat(rand(3,3,1),[1 1 200]);
      r = repmat([1 2 3],[10 1 20]);
      s = pumpkyn.util.multiplyDCM(DCM,r,2);
     sc = (DCM(:,:,1)*r(1,:,1)')';
     ds = s-sc;
     return;
end
           [r,fSeq] = pumpkyn.util.fDim(r,dim3);
                  s = permute(sum(permute(r,[3 2 1]).*DCM,2),[3 1 2]);   
fSeq.postShift(end) = size(s,2);
                  s = pumpkyn.util.eDim(s,fSeq);
end