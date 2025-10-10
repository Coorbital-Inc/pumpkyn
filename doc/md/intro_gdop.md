
# Computing GDOP at the Lunar South Pole

This tutorial demonstrates how to compute the **\*Geometric Dilution of Precision (GDOP)**\* for a surface observer located at the lunar south pole, using a constellation of satellites in a tulip\-shaped orbit family modeled in the Earth–Moon CR3BP.


GDOP quantifies how the geometry of the visible satellites affects position and timing accuracy for navigation or communication users. Lower GDOP values indicate better geometric diversity and higher precision.

# Define an Optimized Satellite Constellation
```matlab
% These time offsets (dtau) distribute satellites throughout a 
% 6:5 resonant, 7-petal tulip-shaped orbit family.
% The configuration was obtained from a global optimization method
% (e.g., genetic algorithm) that minimizes GDOP and maximizes coverage.
dtau =  [0.171755901906622
         0.538936510214872
1.91530114835572
         2.28537063968867
         3.66203724308906
         4.02638956895957]';
   Nr = 1;
   pm = -1;
   Np = 7;
 tau0 = (5/6)*2*pi;
% Generate the constellation using the CR3BP propagation routine.
% Outputs:
%   tau   – Dimensionless time vector
%   rTgt  – Satellite position history [x,y,z] (ND)
%   vTgt  – Satellite velocity history [xdot,ydot,zdot] (ND)
%   mu    – Earth–Moon mass ratio
%   lStar – Characteristic length scale [km]
%   tStar – Characteristic time scale [s]
[tau,rTgt,vTgt,mu,lStar,tStar] = pumpkyn.cr3bp.tulipConstellation(Np, ...
                                 tau0,Nr,pm,dtau);
```
# Interpolate Trajectories based on discretized time step:

Define a fixed time step for interpolation (30 minutes)

```matlab
                      dt = 30*60./tStar;  %timestep [ND]
% Generate an evenly spaced dimensionless time vector for uniform sampling
                   tau_i = linspace(0,tau(end),round(tau(end)/dt))';
% Interpolate both position and velocity vectors using spline interpolation
                    rTgt = interp1(tau,rTgt,tau_i,'spline');
                    vTgt = interp1(tau,vTgt,tau_i,'spline');
                     tau = tau_i;

% This ensures that all satellites are sampled synchronously in time,
% which is essential for GDOP and visibility calculations.
```
# Determine Satellite Visibility from the South Pole Observer
```matlab
               minElAng = 0*pi/180;
                   rObs = [1-mu,0,-1738.1/lStar]; %South Pole
                    rP2 = [1-mu,0,0];
% Compute elevation angles of each satellite relative to the observer.
% phi: elevation angle matrix [time x 1 x Nsats]
                    phi = pumpkyn.cr3bp.elAng(rObs,rTgt,rP2,2);
% Generate a logical visibility mask for satellites below elevation threshold.
% elMask is used to exclude satellites not in view for GDOP computation.
                 elMask = permute(phi < minElAng,[1 3 2]);
                    
```
# Compute the Number of Satellites in View Over Time
```matlab
% Count the number of satellites with elevation above the mask angle
Nsats = sum(phi >= minElAng,3);

% Plot the number of satellites in view as a function of time
figure('Color',[1 1 1]);
stairs(tau.*tStar./86400,Nsats);
set(gca,'YTick', 0:1:numel(dtau)+1); ylim([0 numel(dtau)+1]);
grid on; xlabel('Time [Days]'); ylabel('Satellites In View');
title('Lunar South Pole Visibility from Tulip Constellation');
```
# Interpretation:

Steps indicate discrete changes in which satellites are visible. At least four satellites are needed for 3D navigation. Periodic gaps reveal coverage limits inherent to the resonance geometry.

# Compute the GDOP of an Observer:

Compute the full Dilution of Precision (DOP) metrics. pumpkyn.cr3bp.dop() returns a matrix with GDOP, PDOP, HDOP, and VDOP values.

```matlab
 dop = pumpkyn.cr3bp.dop(rObs,rTgt,elMask,2);
gdop = dop(:,1); % Extract GDOP (overall geometric precision metric)
figure('color',[1 1 1]);
plot(tau.*tStar./86400,gdop);
grid on; ylim([0 8]); xlabel('Time [Days]'); ylabel('GDOP'); 
title('Geometric Dilution of Precision (GDOP) at Lunar South Pole');
```
# Interpretation:

Lower GDOP → better satellite geometry and stronger navigation accuracy. Higher GDOP → satellites are clustered or low on the horizon. The periodic GDOP oscillations reflect the constellation's resonant motion.


This analysis can help identify optimal orbital geometries and phasing for continuous PNT coverage of the lunar south pole.

