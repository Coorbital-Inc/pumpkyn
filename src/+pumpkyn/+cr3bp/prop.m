function [tau,x] = prop(tau,x0,mu,dim3)
%% Purpose:
%
%  This routine will take the initial states of a CR3BP and propagate them
%  forward using a variable precision ode solver to do so.
%
%  Set the matlab setting back to defaults for the symbolic toolbox
%  via sympref('default') if you cannot view the higher precision numbers
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
%  dim3                     integer             Singleton dimension
%                                               specifier.
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
             dim3 = 1;
              phi = eye(6);
               mu = 1.215058560962404E-2;
               x0 = repmat([9.2294247990951594E-1	
                            0.0000000000000000E+0	
                            2.1592634764719765E-1	
                            4.2453615136710976E-13	
                            1.2786666755130519E-1	
                           -2.3343675990673871E-12
                            phi(:)],[1 3]);
             tau0 = 1.8036821222727220E+0; 
          [tau,x] = pumpkyn.cr3bp.prop(tau0,x0,mu,dim3);
          
          figure('color',[1 1 1]);
          plot3(x(:,1,1),x(:,2,1),x(:,3,1),'-r'); hold on;
          plot3(x(:,1,2),x(:,2,2),x(:,3,2),'.b');
          plot3(x(:,1,3),x(:,2,3),x(:,3,3),'ok');
          axis equal;
          
          return;
end

if ~exist('dim3','var')
    dim3 = 2;
end
             opts = odeset('RelTol',1e-10,'AbsTol',1e-10);
          if isscalar(tau)
             tauVec = [0 tau]; 
          else
             tauVec = tau; 
          end

% Flatten the inputs:
               x0 = pumpkyn.util.fDim(x0,dim3);
               x0 = permute(x0,[2 1]);
            oDims = size(x0);
               x0 = x0(:);
% Call the solver:        
          [tau,x] = ode113(@eomND,tauVec,x0(:),opts,mu,oDims);
% Reshape the outputs:
                x = reshape(x,[numel(tau) oDims(1) oDims(2)]);

end

function  etaDot = eomND(tau,eta,mu,oDims)
    eta = reshape(eta,[oDims(1),oDims(2)]);
 etaDot = pumpkyn.cr3bp.eom(tau,eta,mu,1);
 etaDot = etaDot(:);
end
