function Lpts = lagrangePts(mu)
%% Purpose:
%
%  This routine will compute the five lagrange points for any Circular
%  Restricted Three Body Problem given it's mass ratio mu where
%  mu = m2/(m1+m2) such that m1 > m2 >> m3
%
%% Source Reference:
%   Appendix A: Locating The Lagrange Points in Low-Energy Lunar Trajectory 
%   Design. Jeffrey Parker and Rodney Anderson. JPL. 2013.   
%
% Compare against JPL Periodic Orbits Database:
% ---------------------------------------------
%   Earth-Moon System with  mu =  1.215058560962404E-2;
%   L1:	0.83691513	0.00000000	0.00000000
%   L2:	1.15568217	0.00000000	0.00000000
%   L3:	-1.00506265	0.00000000	0.00000000
%   L4:	0.48784941	0.86602540	0.00000000
%   L5:	0.48784941	-0.86602540	0.00000000
%
%   Sun-Earth System with mu = 3.054200000000000E-6
%   L1:	0.98997092	0.00000000	0.00000000
%   L2:	1.01009044	0.00000000	0.00000000
%   L3:	-1.00000127	0.00000000	0.00000000
%   L4:	0.49999695	0.86602540	0.00000000
%   L5:	0.49999695	-0.86602540	0.00000000
%
%% Inputs:
%
%  mu                   [N x 1]             Mass Ratio of the two primary
%                                           bodies, P1 and P2 where
%                                           m1 > m2 and mu = m2/(m1+m2)
%% Outputs:
%
%  Lpts                 [5 x 3 x N]         Lagrange Point Locations in
%                                           the circular restricted
%                                           three-body frame of reference
%                                           Row 1: L1 [x,y,z]
%                                           Row 2: L2 [x,y,z]
%                                           Row 3: L3 [x,y,z]
%                                           Row 4: L4 [x,y,z]
%                                           Row 5: L5 [x,y,z]
%
%
%% Revision History:
%  Darin C. Koblick                                         (c)  10/20/2025
%  Copyright 2025 Coorbital, Inc.
%% ---------------------- Begin Code Sequence -----------------------------
if nargin == 0
    %% Test Solution Accuracy for Earth-Moon and Sun-Earth Systems:
     mu =  permute([1.215058560962404E-2,3.054200000000000E-6],[1 3 2]);
   Lpts = pumpkyn.cr3bp.lagrangePts(mu);
   %Verify Solutions:
       x = Lpts(:,1,:);
       y = Lpts(:,2,:);
       z = Lpts(:,3,:);
    ydot = z;
    xdot = z;
       d = sqrt((x+mu).^2 + y.^2 + z.^2);   
       r = sqrt((x-1+mu).^2 + y.^2 + z.^2);
   xddot = -(1-mu).*(x+mu)./d.^3 - mu.*(x-1+mu)./r.^3 + 2.*ydot+x;
   yddot = -(1-mu).*y./d.^3 - mu.*y./r.^3 - 2.*xdot + y;
   zddot = -(1-mu).*z./d.^3 - mu.*z./r.^3;
   for tm=1:numel(mu)
   fprintf(1,'%s = %0.2E for mu* = %0.2E\n','Solution Accuracy', ...
           max(sqrt(xddot(:,1,tm).^2+yddot(:,1,tm).^2+zddot(:,1,tm).^2)),mu(tm));
   end
   %% Determine the lagrange Point Positions vs various mass ratios:
      muCoef = (-10:0.01:log10(0.95))';
          mu = 10.^muCoef;
    Lpts = pumpkyn.cr3bp.lagrangePts(mu);
    figure('color',[1 1 1]);
    titleStr = {'L_1','L_2','L_3','L_4/L_5'};
    for tl=1:4
        subplot(2,2,tl);
        semilogx(mu,squeeze(Lpts(tl,1,:)),'k','linewidth',3); 
        grid on;
        xlabel('\mu^*');
        ylabel('X position');
        title(titleStr{tl});
    end
   Lpts = [];
   return;
end
%% Set the convergence tolerance:
       tol = 1E-14;
%% Initialize All Points:
      Lpts = zeros(5,3,numel(mu));
%% Reshape Mu to [1 x 1 x N]
        mu = permute(mu(:),[3 2 1]);
%% Compute L4 and L5:
  Lpts(4,1,:) = 0.5-mu;
  Lpts(5,1,:) = Lpts(4,1,:);
  Lpts(4,2,:) = sqrt(3)/2;
  Lpts(5,2,:) = -Lpts(4,2);
%% Initialize all gamma values:
        gamma = NaN(size(mu));
%% Initial Guess for all Gamma values:  
  gammaInit = (mu.*(1-mu)./3).^(1/3);
%% Compute L3:
    gamma0 = gammaInit;
       idx = true(size(gamma0));
while any(idx(:))
    gamma(idx) = ((1-mu(idx)).*(gamma0(idx)+1).^2./ ...
                  (1+2.*mu(idx)+gamma0(idx).*(2+mu(idx)+gamma0(idx)))).^(1/3);
           idx = abs(gamma-gamma0) > tol;
   gamma0(idx) = gamma(idx);
end
   Lpts(3,1,:) = -mu-gamma;  
%% Compute L2:
gamma0 = gammaInit;
   idx = true(size(gamma0));
while any(idx(:))
   gamma(idx) = (mu(idx).*(gamma0(idx)+1).^2./ ...
                (3-2.*mu(idx)+gamma0(idx).*(3-mu(idx)+gamma0(idx)))).^(1/3);
          idx = abs(gamma-gamma0) > tol;
  gamma0(idx) = gamma(idx);
end
  Lpts(2,1,:) = 1-mu+gamma;
%% Compute L1:
gamma0 = gammaInit;
   idx = true(size(gamma0));
while any(idx(:))
   gamma(idx) = (mu(idx).*(gamma0(idx)-1).^2./ ...
                (3-2.*mu(idx)-gamma0(idx).*(3-mu(idx)-gamma0(idx)))).^(1/3);
          idx = abs(gamma-gamma0) > tol;
  gamma0(idx) = gamma(idx);
end
Lpts(1,1,:) = 1-mu-gamma;
end