function [x0,tau0,converged,nullDF,iter,convErr] = cont(x0,tau0,ds,nullDF,mu,tol)
%% Purpose:
%
%  Given an initial state vector in the CR3BP system, this routine
%  will compute a numerical solution such that it has a repeating
%  orbit using psuedo-arclength continuation.
%
%% Inputs:
%
%  x0                       [1 x 6]             Initial Dimensionless
%                                               State Vector of orbit to
%                                               continue
%
%  tau0                     double              Corresponding dimensionless
%                                               Period
%
%
%  ds                       double              step size (+/-)
%
%
%  nullDF                  [7 x 1]              Previous null vector
%                                               gradient (provides
%                                               directional information
%                                               combined with sign of ds)
%
%  mu                       double              Dimensionless Mass Ratio
%                                               of the two primary bodies
%                                               mu = mu2/(mu1+mu2) where
%                                               m1 > m2
%
% tol                       double              Maximum absolute tolerance 
%                                               after propagating a full
%                                               period: norm(H) <tol 
%
% fixPeriod                 integer             true  = tau0 is fixed
%                                               false = tau0 is free
%
%% Outputs:
%
%  x0                       [1 x 6]            Refined Dimensionless 
%                                              State Vector that has a 
%                                              solution such that  
%                                              |x(0)-x(tau0)| < tol
%
% tau0                      double             Refined Dimensionless Period
%
%
% converged                 boolean             true = converged to sol
%                                              false = did not converge
%
% nullDF                    [7 x 1]            Null vector gradient
%
% iter                      integer            Number of iterations req
%
% convErr                   double             norm(H)
%
%
%% Revision History:
%  Darin C. Koblick                                              08/26/2025
%  Copyright 2025 Coorbital, Inc.
%% --------------------------- Begin Code Sequence ------------------------
if nargin == 0
              tol = 1e-12;
               mu = 1.215058560962404E-2;
               x0 = [9.2294247990951594E-1	
                     0.0000000000000000E+0	
                     2.1592634764719765E-1	
                     4.2453615136710976E-13	
                     1.2786666755130519E-1	
                     -2.3343675990673871E-12];
                           ds = 1e-2;
                       nullDF = zeros(7,1);
                         tau0 = 1.8036821222727220E+0;
   [x0,tau0,converged,nullDF,iter] = pumpkyn.cr3bp.cont(x0,tau0,0,nullDF,mu,tol);
                             [~,x] = pumpkyn.cr3bp.prop(tau0,x0,mu);
          figure('color',[1 1 1]);
          plot3(1-mu,0,0,'.k','markersize',10); hold on;
          plot3(x(:,1),x(:,2),x(:,3),'r'); axis equal;
          for tt=1:20
     [x0,tau0,converged,nullDF,iter,convErr] = pumpkyn.cr3bp.cont(x0,tau0, ...
                                               ds,nullDF,mu,tol); 
                        [~,x] = pumpkyn.cr3bp.prop(tau0,x0,mu);
                        plot3(x(:,1),x(:,2),x(:,3),'r');
                        drawnow;
                        fprintf('Number of Iter = %d, Error = %e\n', ...
                                iter,convErr);
          end
          axis equal;
          grid on;
        return;
end
     convErr = 1;
     maxIter = 100;
        iter = 0;
        PHI0 = eye(6);
  nullDFstar = nullDF;
   converged = false;
          x0 = x0(:);
      [~,xx] = pumpkyn.cr3bp.prop([0 tau0],[x0(1:6); PHI0(:)],mu);
          DF = computeFDeriv([x0(1:6); tau0],xx(end,:),mu);
      nullDF = null(DF, tol);
      nullDF = nullDF(:, 1); 
          Vd = [x0(1:6); tau0];
           V = Vd + ds*nullDF;      % Initial guess of the next orbit
           
     if dot(nullDFstar,nullDF,1) < 0
        nullDF = -nullDF; 
     end
     
     while (convErr > tol) && iter < maxIter 
         [~,xx] = pumpkyn.cr3bp.prop([0 V(end)],[V(1:6); PHI0(:)],mu);           
             DF = computeFDeriv(V,xx(end,:),mu);      
             DH = [DF; nullDF(:)'];
              H = updateConstraintVec(xx(1,1:6),xx(end,1:6),V,Vd,nullDF,ds);
        convErr = norm(H);
        
        if convErr > 1e3
            iter = maxIter;
            break;
        end
         V = V - pinv(DH)*H;
      iter = iter + 1;
     end
     
    if (iter < maxIter)
        converged = true;
    end    
    
   [x0,tau0] = deal(V(1:6),V(7));
end

function DFV = computeFDeriv(V,x,mu)
%% Purpose:
%
%  Take the derivative of the F(V) vector
%
       PHI = reshape(x(7:42),6,6);
      xDot = pumpkyn.cr3bp.eom(V(end),x(1:6),mu,2);
       DFV = [PHI(1:4,1:6)-eye(4,6),xDot(1:4)'
              PHI(6,1:5), PHI(6,6)-1.0, xDot(6)
              0,1.0,0,0,0,0,0];
end


function H = updateConstraintVec(X0,Xf,V,Vd,nullDF,ds)
%% Purpose:
%
%  Construct a new augmented constraint vector to limit distance between
%  subsequent solutions via projection onto tangent space:
%
%  H(V) = |  x_f - x_0       |
%         |  y_f - y_0       |
%         |  z_f - z_0       |
%         | dx_f - dx_0      |
%         | dz_f - dz_0      |
%         |     y0           |
%         | (V-V*)^T nu* - ds|
%
      H = zeros(7,1);
 H(1:5) = Xf([1:4,6])- X0([1:4,6]);
   H(6) = X0(2);
   H(7) = dot(V-Vd, nullDF) - ds;
end