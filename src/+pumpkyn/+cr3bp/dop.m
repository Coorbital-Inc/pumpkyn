function dop = dop(rRcv,rSat,maskIdx,dim3)
%% Purpose:
%  
%  this routine will compute various dillusion of precision metrics for
%  an observer given N satellites in view using the barycentered rotating
%  frame (cr3bp).
%   
%% Inputs:
%
%  rRcv                      [t x 3 x 1 x N]        N number of receivers
%
%  rSat                      [t x 3 x M]            M number of Satellites                                                   
%
%
%  maskIdx                   [t x N x M]            Do not consider these
%                                                   satellites when
%                                                   computing the dop
%                                                   metrics as they are
%                                                   masked
%
%  dim3                     integer                 Singleton dimension
%                                                   specifier
%
%% Outputs:
%
%  DOP                       [t x 3 x N]            - Geometric dillution of
%                                                   precision
%                                                   - Position dillution of
%                                                   precision
%                                                   - Time dillution of
%                                                   precision
%
%% Revision History:
%  Darin C. Koblick                                             08-27-2025
%  Copyright 2025 Coorbital, Inc.
%% ---------------------- Begin Code Sequence -----------------------------
if nargin == 0
                   Np = 4;
                 tau0 = 9*2*pi/10;
                   Ns = Np;
                   Nr = 1;
                   pm = -1;
                 dtau = linspace(0,tau0-tau0/Ns,Ns);
  [tau,rSat,~,mu,lStar] = pumpkyn.cr3bp.tulipConstellation(Np,tau0,Nr,pm,dtau);                  
                   rP = 1738.1./lStar;
                 rRcv = [1-mu,0,-rP];
              maskIdx = false(size(rSat,1),size(rSat,3));
                  dop = pumpkyn.cr3bp.dop(rRcv,rSat(:,:,:),maskIdx,2);
     figure('color',[1 1 1])
     plot(tau,dop(:,1,:));
     ylabel('GDOP'); xlabel('Time [ND]');
     grid on;
     ylim([0 20]);
     return;
end
  

%% Assemble the line of site of each satellite-target:
   LOS = rRcv-rSat;       %Line of Sight (receiver - satellite)
   rho = pumpkyn.util.vmag(LOS,dim3);  %Range from satellite(s) to receiver
LOS_uv = LOS./rho;        %Make a Unit Vector of the postions [t x 3 x M x N]
     N = size(LOS,4);     %Number of receivers
     t = size(LOS,1);     %Number of time steps
     
%% Create the cofactor matrix (Qx):  
        qXYZT = NaN(t,4,N);
            
for tr=1:t   %For each time step
  for tt=1:N %For each receiver
           thisLOS = LOS_uv(tr,:,~maskIdx(tr,:,tt),tt);
           if size(thisLOS,3) < 4
               continue;
           end
                 A = ones(size(thisLOS,3),4);
          A(:,1:3) = permute(thisLOS,[3 2 1]);
          %A'*A will alway be a 4 x 4 matrix
          %qXYZT(tr,:,tt) = diag(inv(A'*A));
    qXYZT(tr,:,tt) = inv44Diag(A'*A);
  end
end

dop = cat(dim3,sqrt(sum(qXYZT,2)), ...
               sqrt(sum(qXYZT(:,1:3,:),2)), ...
               sqrt(qXYZT(:,4,:)));

end

function AID = inv44Diag(A)
%% Purpose:
%
%  This routine will compute the inverse of matrix A, which is expected
%  to be a 4 x 4 square matrix and return only the diagonal elements.
%

  a11 = A(1,1,:);  a12 = A(1,2,:); a13 = A(1,3,:); a14 = A(1,4,:);
  a21 = A(2,1,:);  a22 = A(2,2,:); a23 = A(2,3,:); a24 = A(2,4,:);
  a31 = A(3,1,:);  a32 = A(3,2,:); a33 = A(3,3,:); a34 = A(3,4,:);
  a41 = A(4,1,:);  a42 = A(4,2,:); a43 = A(4,3,:); a44 = A(4,4,:);
denom = a11*a22*a33*a44 - a11*a22*a34*a43 - a11*a23*a32*a44 + ...
        a11*a23*a34*a42 + a11*a24*a32*a43 - a11*a24*a33*a42 - ...
        a12*a21*a33*a44 + a12*a21*a34*a43 + a12*a23*a31*a44 - ...
        a12*a23*a34*a41 - a12*a24*a31*a43 + a12*a24*a33*a41 + ...
        a13*a21*a32*a44 - a13*a21*a34*a42 - a13*a22*a31*a44 + ...
        a13*a22*a34*a41 + a13*a24*a31*a42 - a13*a24*a32*a41 - ...
        a14*a21*a32*a43 + a14*a21*a33*a42 + a14*a22*a31*a43 - ...
        a14*a22*a33*a41 - a14*a23*a31*a42 + a14*a23*a32*a41;
term11 =   a22*a33*a44 - a22*a34*a43 - a23*a32*a44 + a23*a34*a42 + a24*a32*a43 - a24*a33*a42;
term22 =   a11*a33*a44 - a11*a34*a43 - a13*a31*a44 + a13*a34*a41 + a14*a31*a43 - a14*a33*a41;
term33 =   a11*a22*a44 - a11*a24*a42 - a12*a21*a44 + a12*a24*a41 + a14*a21*a42 - a14*a22*a41;
term44 =  a11*a22*a33 - a11*a23*a32 - a12*a21*a33 + a12*a23*a31 + a13*a21*a32 - a13*a22*a31;
   AID = abs([term11,term22,term33,term44]./denom);
end
