
# Finding Unstable Manifolds

This tutorial demonstrates how to compute and visualize the **\*invariant manifolds**\* associated with a tulip\-shaped orbit in the Earth–Moon CR3BP.


Invariant manifolds represent natural dynamical pathways: • The **\*stable manifold** **guides trajectories that asymptotically approach the orbit. • The \*\*unstable manifold** \* describes paths that naturally depart from it.


These structures are central to low\-energy transfer design, lunar capture, and disposal trajectory studies. Manifolds also provide insight into orbit stability and long\-term dynamical behavior.


Here, we focus on the **\*unstable manifold**\*, which emerges from small perturbations along the orbit's unstable eigenvectors. Note: orbits with a stability index close to 1 may not generate usable (distinct) manifolds.

# Form the Basis of a Tulip\-Shaped Orbit
```matlab
                     Np = 13;               
                   tau0 = 2*pi; %(5/6)*2*pi;
                    pm = -1; 
% Generate a representative tulip-shaped orbit and its trajectory data.
% Outputs:
%   r      – Position history of the orbit [x,y,z]
%   v      – Velocity history of the orbit [xdot,ydot,zdot]
%   muStar – Earth–Moon mass ratio
%   lStar  – Characteristic length scale [km]
%
% This call returns the full trajectory of a 13-petal periodic orbit.
  [~,r,v,muStar,lStar] = pumpkyn.cr3bp.tulipConstellation(Np,tau0,1,pm,0);
% Select the initial state vector (position + velocity) at t = 0
                    x0 = [r(1,:),v(1,:)];

% Compute the Stability Index
% The stability index k quantifies how perturbations evolve near the orbit.
%   k ≈ 1   → Orbit is neutrally stable
%   k > 1   → Orbit is unstable (has at least one diverging manifold)
%   k < 1   → Orbit is linearly stable
%
% Here, we compute it from the monodromy matrix over one orbital period.
                     k = pumpkyn.cr3bp.stabilityIndex(x0,tau0,muStar);
                     fprintf(1,'Stability Index = %f\n',k);
```
# Determine and Visualize Manifolds

Number of sample points along the orbit which manifolds will be launched

```matlab
        Npts = 6;
      tauVec = linspace(0,tau0-tau0/Npts,Npts);
   [~,x0Vec] = pumpkyn.cr3bp.prop(tauVec,x0,muStar);
% Define the perturbation magnitude (epsilon) for manifold generation.
% The perturbation is scaled inversely by the local velocity magnitude,
% ensuring consistent spatial displacement across the orbit.
     epsilon = (100/(lStar))./sqrt(x0Vec(:,4).^2 + x0Vec(:,5).^2 + x0Vec(:,6).^2);
        
 for tp=1:Npts
     if tp == 1
         pumpkyn.cr3bp.showMoon(lStar,muStar);
         plot3(r(:,1),r(:,2),r(:,3),'-','LineWidth',3,'color',[1 0.5 0]);
         axis equal;
     end
    % Compute the stable (x0_s) and unstable (x0_u) initial conditions
    % by perturbing the orbit along the eigenvectors of the monodromy matrix.
    % The magnitude of the perturbation is set by epsilon(tp).
    [x0_s,x0_u] = pumpkyn.cr3bp.manifolds(tau0,x0Vec(tp,:),muStar,epsilon(tp));
    % Propagate the unstable manifold trajectories forward in time.
    % Each trajectory shows how an infinitesimal perturbation diverges
    % from the periodic orbit naturally under CR3BP dynamics.
    for tm=2
      [tau,x] = pumpkyn.cr3bp.prop(tau0*8,x0_u(tm,:),muStar);
       plot3(x(:,1),x(:,2),x(:,3),'-','Color',[1 1 1 0.3]);
       drawnow;
   end
 end
 plot3(x0Vec(:,1),x0Vec(:,2),x0Vec(:,3),'.','markersize',12,'color','white');
 set(gca,'Clipping','off');
```
# Interpretation:

The white trajectories are the unstable manifolds—natural escape paths. Each originates from a slightly perturbed state along the tulip orbit. Their geometry reveals low\-energy transfer routes and long\-term dynamical structure within the cislunar region.

