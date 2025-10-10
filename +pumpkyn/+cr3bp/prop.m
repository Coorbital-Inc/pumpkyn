function [tau,x] = prop(tau,x0,mu)
%% Purpose:
%
%  This routine will take the initial states of a CR3BP and propagate them
%  forward using a variable precision ode solver to do so.
%
%  set the matlab setting back to defaults for the symbolic toolbox
%  via sympref('default') if you cannot view the higher preicison numbers
%  in your command window properly.
%
%% Inputs:
%
%  tau                      [N x 1]             Dimensionless Time Vector
%                                               to propagate the initial
%                                               state, x0, corresponding to
%                                               tau0 = 0
%
%  x0                       [1 x 42]            Dimensionless state
%                                               corresponding to tau0 = 0
%
%  mu                       double              mass ratio of the two
%                                               primaries mu =
%                                               mu2/(mu1+mu2)
%
%% Outputs:
%
%  tau                      [N x 1]               Dimensionless Time Vector
%
%   x                       [N x 42]              Dimensionles State Vector
%                                                 [x,y,z,xdot,ydot,zdot,PHI]
%
%% Revision History:
%  Darin C. Koblick                                              08/26/2025
%  Copyright 2025 Coorbital, Inc.
%% --------------------------- Begin Code Sequence ------------------------



if nargin == 0
              phi = eye(6);
               mu = 1.215058560962404E-2;
               x0 = [9.2294247990951594E-1	
                     0.0000000000000000E+0	
                     2.1592634764719765E-1	
                     4.2453615136710976E-13	
                     1.2786666755130519E-1	
                     -2.3343675990673871E-12
                     phi(:)];
             tau0 = 1.8036821222727220E+0; 
          [tau,x] = pumpkyn.cr3bp.prop(tau0,x0,mu);
          
          figure('color',[1 1 1]);
          plot3(x(:,1),x(:,2),x(:,3),'-r'); hold on;
          axis equal;
          
          return;
end

             opts = odeset('RelTol',1e-10,'AbsTol',1e-10);
          if isscalar(tau)
             tauVec = [0 tau]; 
          else
             tauVec = tau; 
          end
          [tau,x] = ode113(@pumpkyn.cr3bp.eom,tauVec,x0,opts,mu,1);         
end