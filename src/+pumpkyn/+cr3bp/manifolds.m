function [x0_s,x0_u] = manifolds(tau0,x0,muStar,epsilon)
%% Purpose:
%
%  Given an inital state of an unstable peridic three-body orbit, this 
%  routine will compute both the stable and unstable invariant manifold 
%  initial states. This will be based off a specified perturbation defined
%  by a dimensionless length
%   
%  epsilon = length / sqrt(vx^2 + vy^2 + vz^2)
%
%
%% Inputs:
%
%  tau0                 double              Dimensionless Period of Orbit
%
%  x0                   [1 x 6]             Dimensionless States (pos/vel)
%
%  muStar               double              Mass Ratio of Primaries
%                                           muStar = mu2/(mu1+mu2)
%
%  epsilon              double              perturbation parameter
%                                           
%
%
%% Ouputs:
%
%  x0_s                 [N x 6]             Initial dimensionless states of 
%                                           the Stable Invariant Manifolds
%                                           Includes Interior and exterior
%                                           directions
%
%  x0_u                 [N x 6]             Initial dimensionless states of 
%                                           the Unstable Invariant Manifolds
%                                           Includes Interior and exterior
%                                           directions
%                                            
%% Revision History:
%  Darin C. Koblick                                              09/15/2025
%  Copyright 2025 Coorbital, Inc.
%% --------------------------- Begin Code Sequence ------------------------
if nargin == 0
   %two-dimensional lyaponuv orbit bout L2 
   tau0 = 3.4344774781239868;
     x0 = [1.1130439425606091E+0	
           -6.0030029039812299E-28	
           8.9831983586449705E-33	
           4.0854090008117549E-15	
           2.1082532553834954E-1	
           7.5485662781839097E-33]';
   muStar = 1.215058560962404E-2;
    tStar = 382981;
    lStar = 389703;
  
  %Propagate orbit to each point:
      Npts = 90;
    tauVec = linspace(0,tau0-tau0/Npts,Npts);
 [~,x0Vec] = pumpkyn.cr3bp.prop(tauVec,x0,muStar);
    figure('color',[1 1 1]);
    plot3(x0Vec(:,1),x0Vec(:,2),x0Vec(:,3),'.k','markersize',14); hold on;     
  epsilon = (100/(lStar))./sqrt(x0Vec(:,4).^2 + x0Vec(:,5).^2 + x0Vec(:,6).^2);
    
  
  %Propagate manifolds for each point
  for tp = 1:Npts
  [x0_s,x0_u] = pumpkyn.cr3bp.manifolds(tau0,x0Vec(tp,:),muStar,epsilon(tp));
  for tm=1:2
      [tau,x] = pumpkyn.cr3bp.prop(tau0*1.5,x0_u(tm,:),muStar);
       plot3(x(:,1),x(:,2),x(:,3),'r');
  end
  end
   axis equal;
   set(gca,'Clipping','off');

  return;
end

%% Get the Stbale and Unstable Eigen Vectors:
        data = pumpkyn.cr3bp.orbitProperties(x0,tau0,muStar);
         V_s = data.StableEigenVecs;
         V_u = data.UnstableEigenVecs;

%% Initial Conditions for stable and unstable manifolds:
         eps_V_s = real(epsilon.*V_s./pumpkyn.util.vmag(V_s,2));
         eps_V_u = real(epsilon.*V_u./pumpkyn.util.vmag(V_u,2));
            x0_s = [x0 + eps_V_s; x0 - eps_V_s];  %Interior; Exterior        
            x0_u = [x0 + eps_V_u; x0 - eps_V_u];  %Interior; Exterior
end