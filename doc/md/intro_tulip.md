
# Creating a Tulip\-Shaped Orbit

This tutorial demonstrates how to: 1. Generate a family of tulip\-shaped periodic orbits in the Earth–Moon CR3BP 2. Compute and visualize key dynamical properties for each orbit 3. Understand how these properties evolve with orbital period


Tulip\-shaped orbits are resonant three\-body trajectories characterized by repeating petal\-like patterns around the Moon. Each family is defined by its number of petals (Np) and hemisphere (pm), representing distinct energy and symmetry configurations.

# Retrieve All Tulip\-Shaped Orbits for a Specified Petal Count
```matlab
                       rMoon = 1737.1; % Lunar radius [km]
                          Np = 7;      % Number of petals
                          pm = +1;     % +1 = Northern, -1 = Southern
% Retrieve all cataloged orbits for the selected family.
% The function returns:
%  tau0  – Dimensionless orbital period of each orbit
%  x0    – Corresponding initial state vectors [x, y, z, xdot, ydot, zdot]
%  mu    – Earth–Moon mass ratio
%  tStar – Characteristic time scaling (s)
%  lStar – Characteristic length scaling (km)
[tau0, x0, mu, tStar, lStar] = pumpkyn.cr3bp.getTulip([],Np,pm);
```
# Compute Properties for Each Family Member

Preallocate storage arrays for performance

```matlab
StabilityIndex = zeros(size(x0, 1), 1);
        Jacobi = zeros(size(x0, 1), 1);
      Perilune = zeros(size(x0, 1), 1);
       Apolune = zeros(size(x0, 1), 1);
  DoublingTime = zeros(size(x0,1),  1);
   maxLunarOcc = zeros(size(x0,1),  1);
   totLunarOcc = zeros(size(x0,1),  1);

for ts=1:size(x0,1)
                    fprintf(1,'Processing %d of %d\n',ts,size(x0,1));
                    % Compute orbit properties
                    data = pumpkyn.cr3bp.orbitProperties(x0(ts,:),tau0(ts),mu,lStar);
    StabilityIndex(ts,1) = data.StabilityIndex;
            Jacobi(ts,1) = data.Jacobi;
          Perilune(ts,1) = data.Perilune.*lStar - rMoon;
           Apolune(ts,1) = data.Apolune.*lStar - rMoon;
      DoublingTime(ts,1) = data.DoublingTime.*tStar./86400;
       maxLunarOcc(ts,1) = data.MaxLunarOcc.*tStar./3600;
       totLunarOcc(ts,1) = data.TotLunarOcc.*tStar./3600;
end
```
# Plot Doubling Time vs. Period

Shows how long it takes for small perturbations to double in magnitude. Longer doubling times generally indicate greater orbital stability.

```matlab
figure('Color',[1 1 1]);
plot(tau0.*tStar./86400,DoublingTime);
grid on;
xlabel('Period [Days]'); ylabel('Doubling Time [Days]');
```
# Plot Stability Index vs. Period

The stability index measures how sensitive an orbit is to perturbations. Values near 1 imply neutral stability, while larger values indicate increasingly unstable behavior.

```matlab
figure('Color',[1 1 1]);
plot(tau0.*tStar./86400,StabilityIndex);
grid on;
xlabel('Period [Days]'); ylabel('Stability Index');
```
# Plot Perilune and Apolune Altitude vs. Period

Visualizes how the minimum and maximum lunar altitudes change as a function of orbital period. This helps identify orbits with practical altitudes for communications or observation missions.

```matlab
figure('Color',[1 1 1]);
plot(tau0.*tStar./86400,Perilune./1000,'r'); hold on;
plot(tau0.*tStar./86400,Apolune./1000,'b');
grid on;
xlabel('Period [Days]'); ylabel('Surface Alt [Mm]');
legend('Perilune','Apolune');
```
# Plot Jacobi Constant vs. Period

The Jacobi constant represents the total energy of each orbit in the rotating frame. It provides insight into the energy distribution and dynamical accessibility of each family member.

```matlab
figure('Color',[1 1 1]);
plot(tau0.*tStar./86400,Jacobi);
grid on;
xlabel('Period [Days]'); ylabel('Jacobi Constant');
```
# Plot Max Lunar Occultation vs. Period

Displays how long each orbit spends behind the Moon relative to Earth, representing communications blackout duration. Shorter occultations are preferred for continuous Earth visibility.

```matlab
figure('Color',[1 1 1]);
plot(tau0.*tStar./86400,maxLunarOcc);
grid on;
xlabel('Period [Days]'); ylabel('Max Lunar Occultation [Hrs]');
```
