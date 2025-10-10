%% Purpose:
%
% A zero-velocity surface is a concept that relates to the N-body problem 
% of gravity. It represents a surface a body of given energy cannot cross, 
% since it would have zero velocity on the surface.

%% Compute the Jacobi Constant for all points in Cislunar Space
muStar = 1.215058560962404E-2;
 [x,y] = meshgrid(0.5:0.001:1.2,-0.3:0.001:0.3);
 zeros = x.*0;
    x0 = [x(:),y(:),zeros(:),zeros(:),zeros(:),zeros(:)];
     J = pumpkyn.cr3bp.jacobi(x0,muStar);
     J = reshape(J,size(x));
J(J>=4) = NaN;

%%  Show the contour Plot:
figure('Color',[1 1 1]);
contourf(x,y,J);
set(gca,'ColorScale','log'); hold on;
c = colorbar;
c.Label.String = 'Jacobi Constant';
xlabel('x [ND]'); ylabel('y [ND]');


%% Interpretation of Zero-Velocity Surface (ZVS)
%
% The Jacobi constant for this 9-petal tulip-shaped orbit is C = 3.1777.
% In the Earth–Moon CR3BP, this value lies BETWEEN the L1 and L2 thresholds
% (C_L1 ≈ 3.188, C_L2 ≈ 3.172).
%
% → Because C < C_L1, the neck near L1 is OPEN, allowing potential
%   transfer trajectories from the lunar region toward Earth.
%
% → However, since C > C_L2, the neck near L2 remains CLOSED,
%   preventing natural escape into the exterior (deep-space) region.
%
% In summary, this orbit exists in a transitional energy regime where
% motion is energetically permitted through L1 toward Earth, but still
% confined against escape through L2 into interplanetary space.
                        tau0 = 2*pi;
                          Np = 9;
                          pm = +1;
                  [tau0, x0] = pumpkyn.cr3bp.getTulip(tau0,Np,pm);
                    [tau, x] = pumpkyn.cr3bp.prop(tau0,x0,muStar);
                          J0 = pumpkyn.cr3bp.jacobi(x0,muStar);
plot(x(:,1),x(:,2),'r','linewidth',2);
title(['\color{red}Jacobi Constant = ',num2str(J0)]);
axis equal; grid on;                          

                   