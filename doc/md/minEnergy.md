
# Minimum Energy Transfer to the Lunar Surface

This demonstration shows how to use Pumpkyn's CR3BP routines to estimate the minimum energy (ΔV) required for a spacecraft in a tulip\-shaped orbit to reach the lunar surface—or equivalently, the energy required for ascent from the surface to that orbit.


The method evaluates instantaneous transfer opportunities between points along a periodic orbit and a grid of surface locations, then determines the minimum achievable ΔV across the entire surface.

# Extract an NRHO 9:2 Resonant Orbit (Low\-Energy Example)
```matlab
                            tau0 = (2/9)*2*pi;
                              Np = 1;
                              pm = -1;
[tau0, x0, muStar, tStar, lStar] = pumpkyn.cr3bp.getTulip(tau0,Np,pm);
```
# Propagate the Orbit for One Full Revolution

This provides the full periodic trajectory in the rotating frame, which serves as the reference path for computing transfer opportunities.

```matlab
                         [tau,x] = pumpkyn.cr3bp.prop(tau0, x0, muStar);
```
# Generate a Grid of Points Representing the Lunar Surface
```matlab
                          rP = 1738.1./lStar;
                [xf,lla,xyz] = pumpkyn.cr3bp.pointSphere(500,[1-muStar,0,0],rP);
                          xf = [xf,xf.*0];
% The sphere is discretized into 500 unique surface points centered at the Moon's 
% barycentric position. These points serve as potential landing sites.
```
# Compute the Minimum ΔV to Reach Each Surface Point
```matlab
                       dVmin = NaN(size(xf,1),1);
% For each point along the orbit, compute the required ΔV to reach every 
% surface point. The minimum across all orbital phases gives the best 
% (lowest energy) transfer option to that location.
for tt=1:numel(tau)
      dV = pumpkyn.cr3bp.minDeltaV(repmat(x(tt,:),[size(xf,1) 1]),xf,muStar);
   dVmin = min([dVmin,dV],[],2);
end
% Reshape into 2D surface grid and convert to physical units [km/s]
   dVmin = reshape(dVmin,[size(xyz,1),size(xyz,2)]).*lStar./tStar; %km/s
```
# Plot the Minimum Energy ΔV Distribution (NRHO Case)
```matlab
figure('color',[0 0 0]);
surface(xyz(:,:,1),xyz(:,:,2),xyz(:,:,3),dVmin); hold on;
set(gcf,'color','k'); axis off; hc = colorbar;
plot3(x(:,1),x(:,2),x(:,3),'w','linewidth',2);
axis equal;
shading interp; set(hc,'color','w');
hc.Label.String = '\color{white}Min \DeltaV_{tot} [km/s]'; view(270,0);
title('\color{white}Minimum Energy \DeltaV Lunar Surface Access');
set(gca,'clipping','off');
```
# Extract a Higher\-Energy 15\-Petal Tulip Orbit
```matlab
                            tau0 = 2*pi;
                              Np = 15;
[tau0, x0, muStar, tStar, lStar] = pumpkyn.cr3bp.getTulip(tau0,Np,pm);
```
# Propagate the Orbit for One Full Revolution
```matlab
                         [tau,x] = pumpkyn.cr3bp.prop(tau0, x0, muStar);
```
# Recompute the Minimum ΔV Map for the 15\-Petal Orbit
```matlab
                       dVmin = NaN(size(xf,1),1);
for tt=1:numel(tau)
      dV = pumpkyn.cr3bp.minDeltaV(repmat(x(tt,:),[size(xf,1) 1]),xf,muStar);
   dVmin = min([dVmin,dV],[],2);
end
   dVmin = reshape(dVmin,[size(xyz,1),size(xyz,2)]).*lStar./tStar; %km/s
```
# Plot the Minimum Energy ΔV Distribution (Tulip Orbit)
```matlab
figure('color',[0 0 0]);
surface(xyz(:,:,1),xyz(:,:,2),xyz(:,:,3),dVmin); hold on;
set(gcf,'color','k'); axis off; hc = colorbar;
plot3(x(:,1),x(:,2),x(:,3),'w','linewidth',2);
axis equal;
shading interp; set(hc,'color','w');
hc.Label.String = '\color{white}Min \DeltaV_{tot} [km/s]'; view(270,0);
title('\color{white}Minimum Energy \DeltaV Lunar Surface Access');
set(gca,'clipping','off');
```
# Interpretation:

Color shading indicates the minimum ΔV required to reach each surface point. Low ΔV regions correspond to natural low\-energy access zones. Comparing different orbits (NRHO vs. Tulip) highlights how orbital energy and geometry influence feasible landing or ascent corridors.

