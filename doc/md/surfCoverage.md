
# Surface Coverage Analysis

This example demonstrates how to use Pumpkyn's CR3BP tools to perform preliminary lunar surface coverage analysis. The workflow supports both single satellites and multi\-satellite constellations, provided that all vehicle trajectories share a common epoch.

# Initialize a Tulip Constellation

Generate a constellation of satellites in tulip\-shaped orbits using the cr3bp.tulipConstellation() routine. Parameters include: Np – number of tulip petals p,q – resonance ratio (Moon revolutions vs. Earth revolutions) Nr – number of orbit periods to propagate pm – hemisphere selection (+1 = Northern, –1 = Southern apolune) Ns – number of satellites in constellation dtau – relative offsets in pseudo\-anomaly (dimensionless)

```matlab
                 Np = 15;       
                  p = 1;        
                  q = 1;        
                 Nr = 2;       
                 pm = +1;       
               tau0 = q*2*pi/p;
                 Ns = 3;       
               dtau = linspace(0,0.5*tau0/Ns,Ns); %Pseudo-Anomaly (ND)
 [tau,r,v,muStar,lStar,tStar] = pumpkyn.cr3bp.tulipConstellation(Np,tau0,Nr,pm,dtau);
```
# Generate Lunar Surface Points

Approximate the Moon's surface with a discrete point sphere. These locations serve as ground sites for coverage evaluation.

```matlab
             rP = 1738.1./lStar;  % Dimensionless lunar radius
           Npts = 500;            % Number of surface points
           rPt0 = [1-muStar,0,0]; % Moon barycentric offset
 [rPts,lla,xyz] = pumpkyn.cr3bp.pointSphere(Npts,rPt0,rP);
```
# Compute Maximum Gap Times

For each surface point, compute the maximum continuous time interval during which no satellite exceeds a specified elevation mask angle.

```matlab
       minElAng = 5*pi/180;                   % Minimum elevation angle [rad]
           dTau = 60*60/tStar;                % Time step [ND]
    maxGapTimes = pumpkyn.cr3bp.maxGaps(tau,r,rPts,rPt0,minElAng,1,dTau);
    maxGapTimes = reshape(maxGapTimes,[size(xyz,1), size(xyz,2)]);
```
# Visualize Coverage Results

Plot the Moon in a body\-fixed frame, color\-coded by maximum coverage gap (hours). Satellite trajectories are overlaid for context.

```matlab
figure('color',[0 0 0]);
surface(xyz(:,:,1),xyz(:,:,2),xyz(:,:,3),maxGapTimes.*tStar./3600); hold on;
set(gcf,'color','k'); axis off; hc = colorbar(); colormap(jet);
plot3(r(:,1,1),r(:,2,1),r(:,3,1),'w');
plot3(squeeze(r(1,1,:)),squeeze(r(1,2,:)),squeeze(r(1,3,:)), ...
      '.','MarkerSize',15);
axis equal; shading interp; set(hc,'color','w');
hc.Label.String = '\color{white}Max Coverage Gap [Hrs]'; view(270,0);
set(gca,'clipping','off');
```
