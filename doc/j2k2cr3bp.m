%% Converting a Low Lunar Orbit (LLO) from J2K to CR3BP Coordinates
%
% This tutorial demonstrates how to convert a low lunar orbit defined in the
% Moon-centered inertial (J2K) frame into the barycentric rotating frame
% used in the Circular Restricted Three-Body Problem (CR3BP).
%
% The conversion allows comparison between traditional two-body orbits 
% (as used in standard mission analysis) and their equivalent representation
% in the CR3BP framework, enabling hybrid modeling and low-energy transfer design.

%% Establish two-body and three-body constants
            jd0 = juliandate(2020,1,1);     %Inital Epoch
             mu = 0.012150585609624;        %Mass ratio
              G = 6.67384e-20;              %Gravitational Constant
              M = 5.9736E24 + 7.35E22;      %Char mass (sum of primaries)

% These constants define the Earth–Moon CR3BP system parameters.
% The gravitational influence of the spacecraft is neglected (restricted case).

%% Define a Low Lunar Orbit (LLO) in the Moon-Centered Inertial Frame
          muMoon = G*mu*M;
           rMoon = 1738.1;
              a0 = rMoon+100;
% Define an initial state vector in the J2K frame:
%   Position  [x, y, z] = [a0, 0, 0]
%   Velocity  [vx, vy, vz] = [0, sqrt(muMoon/a0), 0]
             rv0 = [a0,0,0,0,sqrt(muMoon/a0),0];
% Compute orbital period using Kepler's third law for two-body motion.
               P = 2*pi*sqrt((a0^3)/muMoon);
            
% This represents a 100 km circular low lunar orbit, typical of mapping
% or reconnaissance missions (e.g., LRO-type trajectories).

%% Convert from J2K to CR3BP Frame

% The fromJ2K() routine transforms the initial J2K state into the
% dimensionless CR3BP rotating frame.
%
% Inputs:
%   jd0 – Epoch (Julian Date)
%   rv0 – State vector in J2K frame [km, km/s]
%   mu  – Earth–Moon mass ratio
%   M   – Combined Earth–Moon mass [kg]
%   2   – Primary body identifier (2 = Moon-centered)
%
% Outputs:
%   x0    – Dimensionless CR3BP state vector
%   tStar – Characteristic time scaling [s]
      [x0,tStar] = pumpkyn.cr3bp.fromJ2K(jd0,rv0,mu,M,2);
% Compute the orbit period in CR3BP dimensionless time units
            tau0 = P/tStar;

%% Propagate the Orbit in the CR3BP Frame
% Propagate for several orbital periods to observe motion in the rotating frame.
          [tau,x] = pumpkyn.cr3bp.prop(tau0*5,x0,mu);

%% Visualize the Orbit Relative to the Moon
% Display the Moon at its barycentric position [1 - mu, 0, 0].
pumpkyn.cr3bp.showMoon();
% Plot the transformed orbit trajectory in the rotating frame.
plot3(x(:,1),x(:,2),x(:,3),'b');
axis equal; set(gca,'clipping','off');

%% Interpretation:
% The blue trajectory represents the same 100 km LLO, now expressed 
% in the Earth–Moon rotating reference frame. This conversion enables 
% hybrid modeling where two-body and three-body effects can be 
% compared directly.
