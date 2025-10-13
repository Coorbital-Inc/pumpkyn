%% Purpose:
%
% Demonstrates how to generate and visualize a constellation of satellites 
% evenly spaced in time along a tulip-shaped orbit in the Earth–Moon 
% Circular Restricted Three-Body Problem (CR3BP).
%
% Each satellite follows the same tulip-shaped trajectory, but with 
% different phase offsets (pseudo-anomalies) so that they are 
% uniformly distributed around the orbit. This configuration can 
% represent a synchronous cislunar communications or observation 
% constellation.

%% Initialize the Tulip-Shaped Orbit Constellation:
                 Np = 7;        % Number of petals
                  p = 5;        % Revolutions around moon
                  q = 4;        % Revolutions around earth
                 Nr = 2;        % Number of revolutions to propagate
                 pm = -1;       % +1 = Northern, -1 = Southern Apolune
               tau0 = q*2*pi/p; %Orbital Period
                 Ns = Np;       %Number of Satellites
               dtau = (linspace(0,Ns-1,Ns)).*tau0./Np;  %Pseudo-Anomaly (ND)
% Generate the constellation using the CR3BP propagation routine.
% Outputs:
%  - tau: dimensionless time vector
%  - r: position of each satellite [x,y,z] (dimensionless)
%  - mu: mass ratio parameter (Earth–Moon)
%  - lStar: characteristic length scale (Earth–Moon distance, km)
 [tau,r,~,mu,lStar] = pumpkyn.cr3bp.tulipConstellation(Np,tau0,Nr,pm,dtau);
     
%% Show the Moon-Centered Inertial Trajectory:
% Display a 3D Moon sphere located at its barycentric position [1-mu, 0, 0].
% This helps visualize the orbits relative to the Moon in the rotating frame.
pumpkyn.cr3bp.showMoon(lStar,mu);
% Plot the reference tulip-shaped orbit of the first satellite.
plot3(r(:,1,1),r(:,2,1),r(:,3,1),'w','LineWidth',2);

%% Show Each Starting Position of the Satellites:
% Plot the initial positions of all satellites at t = 0 to show
% how they are evenly distributed along the orbit.
plot3(squeeze(r(1,1,:)),squeeze(r(1,2,:)),squeeze(r(1,3,:)), ...
     '.','MarkerSize',20);
set(gca,'clipping','off');