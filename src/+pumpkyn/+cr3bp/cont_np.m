function [x0,iter] = cont_np(x0,tau0,mu,tol)
%% Purpose:
%
%  Fixed Period Natural Parameter Continuation.  Given an initial guess
%  on dimensionless position and velocity states, this routine will use
%  natural parameter continuation to refine them such that the orbit
%  is repeating given a desired period.
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
%  mu                       double              Mass ratio of the two
%                                               primary bodies
%                                               mu = m2/(m1+m2)
%
%  tol                      double              Convergence tolerance
%
%% Outputs:
%
%  x0                       [1 x 6]             Full Dimensionless States
%                                               of the Tulip-shaped orbit
%                                               corresponding to
%                                               [x,y,z,xdot,ydot,zdot]
%
%% Revision History:
%  Darin C. Koblick                                          (c) 09/30/25
%  Copyright 2025 Coorbital, Inc.
%% ------------------------ Begin Code Sequence ---------------------------
if nargin == 0
              tol = 1e-10;
               mu = 1.215058560962404E-2;
               x0 = [9.2294247990951594E-1	
                     0.0000000000000000E+0	
                     2.1592634764719765E-1	
                     4.2453615136710976E-13	
                     1.2786666755130519E-1	
                     -2.3343675990673871E-12];
            tau0 = 1.8036821222727220E+0;
            dTau = 1e-3;
                         [x0,iter] = pumpkyn.cr3bp.cont_np(x0,tau0,mu,tol); 
                             [~,x] = pumpkyn.cr3bp.prop(tau0,x0,mu);
          figure('color',[1 1 1]);
          plot3(1-mu,0,0,'.k','markersize',10); hold on;
          plot3(x(:,1),x(:,2),x(:,3),'r'); axis equal;
          for tt=1:100
                         tau0 = tau0+dTau;
                    [x0,iter] = pumpkyn.cr3bp.cont_np(x0,tau0,mu,tol); 
                        [~,x] = pumpkyn.cr3bp.prop(tau0,x0,mu);
                        plot3(x(:,1),x(:,2),x(:,3),'r');
                        drawnow;
                        fprintf('Number of Iter = %d\n',iter);
          end
          axis equal;
          grid on;
        return;
end


tau_n = tau0/2;
 PHI0 = eye(6);
for iter=1:100
    %Propagate orbit for half period:
    [~,x] = pumpkyn.cr3bp.prop([0 tau_n],[x0(:); PHI0(:)],mu);
      PHI = reshape(x(end,7:42),[6 6]);
       DF =  [PHI(4,1),PHI(4,3),PHI(4,5);  %xDot == 0
              PHI(6,1),PHI(6,3),PHI(6,5);  %zDot == 0
              PHI(2,1),PHI(2,3),PHI(2,5)]; % y == 0
    x_star = [x0(1); x0(3); x0(5)];      
    x_star = x_star - pinv(DF)*[x(end,4); x(end,6); x(end,2)];
        dy = [x_star(1); x0(2); x_star(2); x0(4); x_star(3); x0(6)] - x0(:);
        x0 = x0(:) + dy;
    %Break out if the difference is within the specified tolerance:
    if norm(dy) < tol
        break;
    end
end
%Reshape output:
x0 = x0(:)';
end