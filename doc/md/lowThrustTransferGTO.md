
# Low\-Thrust Minimum\-Time Transfer: GTO → Tulip\-Shaped Orbit

This live script demonstrates a **\*minimum\-time, low\-thrust transfer** **of a 15 kg CubeSat (Tmax = 25 mN) from a geostationary transfer orbit (GTO\-like, Earth\-centered) into a \*\*7\-petal tulip\-shaped orbit** \* in the Earth–Moon Circular Restricted Three\-Body Problem (CR3BP).


The target tulip orbit exhibits **\*zero lunar occultation****, making it well\-suited for \*\*PNT** **(Positioning, Navigation, and Timing) and \*\*lunar communications relay architectures** \*.


\-\-\- **\*Sections**\* 1. Define constants and nondimensional scaling 2. Set initial GTO state and convert to CR3BP coordinates 3. Fetch tulip seed orbit and propagate 4. Select matching departure and arrival points 5. Define low\-thrust spacecraft and solver parameters 6. Solve the minimum\-time transfer problem 7. Propagate and visualize the transfer trajectory 8. Compute propellant and ΔV results


\-\-\- **\*Dependencies:** **All CR3BP routines are part of the `pumpkyn.cr3bp.** ` API. Variables ending in `_ND` are nondimensional (CR3BP\-scaled).

# Characteristic Quantities and Constants

Earth–Moon system constants used for nondimensional scaling.

```matlab
M = 5.9736E24 + 7.35E22;        % Combined mass of primaries (kg)
G = 6.67384e-20;                % Gravitational constant (km^3/kg/s^2)
muStar = 0.012150585609624;     % Mass ratio μ = m₂ / (m₁ + m₂)
lStar  = 389703.264829278;      % Characteristic length (km)
tStar  = 382981.289129055;      % Characteristic time (s)
```
# Define Initial GTO Orbit (Dimensional → CR3BP ND)

Approximate GTO parameters and conversion to nondimensional coordinates.

```matlab
muEarth = G * (1 - muStar) * M;
rEarth  = 6378;                 % Earth radius (km)
rP = rEarth + 350;              % Perigee altitude (km)
rA = rEarth + 35786;            % Apogee altitude (km)
a0 = (rA + rP) / 2;             % Semi-major axis (km)
e0 = (rA - rP) / (rA + rP);     % Eccentricity
i0 = 0;                         % Inclination (rad)
O0 = 0;                         % RAAN (rad)
o0 = -25 * pi / 180;            % Argument of perigee (rad)
nu0 = 0;                        % True anomaly (rad)

oev0 = [a0, e0, i0, o0, O0, nu0];
[r0, v0] = pumpkyn.cr3bp.orb2eci(muEarth, oev0, 2);
rv0 = [r0, v0];

% Keplerian period (s)
P0 = 2 * pi * sqrt(a0^3 / muEarth);

% Convert to CR3BP nondimensional rotating barycentric frame
rv0ND = pumpkyn.cr3bp.fromPCI(0, rv0, muStar, tStar, lStar, 1);
```
# Target Tulip Orbit Seed

Retrieve a 7\-petal southern\-hemisphere tulip orbit seed.

```matlab
Np = 7;               % Number of petals
pm = -1;              % Southern hemisphere
tau0 = (5/6) * 2*pi;  % Nondimensional period estimate
[~, x0] = pumpkyn.cr3bp.getTulip(tau0, Np, pm, 1e-12);
```
# Propagate Orbits to Identify Transfer Geometry

Propagate both GTO and tulip orbits for one period.

```matlab
[tauTgt, rvTgt] = pumpkyn.cr3bp.prop(tau0, x0, muStar);
[tauInt, rvInt] = pumpkyn.cr3bp.prop(P0 / tStar, rv0ND, muStar);

% Choose departure and arrival states heuristically.
[~, idx_f] = max(rvTgt(:,5));   % Arrival index (e.g., maximum z-dot)
idx_0 = 1;                      % Departure index

rv0 = rvInt(idx_0, :);
rvf = rvTgt(idx_f, :);
```
# Define Low\-Thrust Spacecraft and Solver Parameters

Convert thrust and Isp to nondimensional quantities.

```matlab
m0 = 15;                                   % Initial mass (kg)
g0 = 9.80665 * tStar^2 / (1000 * lStar);   % Gravity in ND units
Tmax = 0.025;                              % Max thrust (N)
Tmax = (Tmax / m0) * tStar^2 / (lStar * 1000); % ND thrust accel
Isp = 2100 / tStar;                        % ND specific impulse
c = Isp * g0;                              % ND exhaust velocity
```
# Solve Minimum\-Time Transfer

Use indirect optimal control to find the minimum\-time trajectory.

```matlab
lambda0_tf = [ 190.476497248065
              -79.7064866984696
              -0.430399154713168
              0.301159446575878
              0.586671892449694
            -0.00711582435720301
              4.32931089137559
              6.29081541876621];

sol_lambda0_tf = pumpkyn.cr3bp.tfMin(rv0, rvf, lambda0_tf, Tmax, c, muStar);
```
# Propagate the Optimized Trajectory

Integrate the resulting low\-thrust trajectory.

```matlab
[tau, rv] = pumpkyn.cr3bp.tfMinProp(sol_lambda0_tf(8), ...
              [rv0, 1, sol_lambda0_tf(1:7)'], Tmax, c, muStar);
```
# Propellant Usage and ΔV Analysis

Compute propellant consumption and cumulative ΔV.

```matlab
mTot = rv(:,7) .* m0;                       % Mass vs. time (kg)
mProp = m0 - mTot;                          % Propellant used (kg)
dVtot = c .* log(m0 ./ mTot) .* lStar ./ tStar;  % ΔV (km/s)

figure('Color',[0 0 0]);
subplot(1,2,1);
plot(tau * tStar / 86400, mProp, 'w', 'LineWidth', 2);
xlabel('Transfer Time [days]'); ylabel('Propellant [kg]');
grid on; set(gca,'Color','k','XColor','w','YColor','w');

subplot(1,2,2);
plot(tau * tStar / 86400, dVtot, 'w', 'LineWidth', 2);
xlabel('Transfer Time [days]'); ylabel('\DeltaV_{tot} [km/s]');
grid on; set(gca,'Color','k','XColor','w','YColor','w');
```
# Visualize Earth, Moon, Tulip, and Transfer

Render the 3D geometry of the transfer.

```matlab
[tauP,rvP] = pumpkyn.cr3bp.prop(tau0,rv(end,1:6),muStar);
        rv = [rv(:,1:6); rvP(2:end,:)];
       tau = [tau; tau(end)+tauP(2:end)];

hIn = figure('color', [0 0 0]);
pumpkyn.cr3bp.showEarth(lStar, muStar, hIn);
pumpkyn.cr3bp.showMoon(lStar, muStar, hIn);  
set(gca, 'color', 'k'); hold on;

Lpts = pumpkyn.cr3bp.lagrangePts(muStar);
plot3(Lpts(1,1), Lpts(1,2), Lpts(1,3), '.b', 'MarkerSize', 20);
text(Lpts(1,1), Lpts(1,2), Lpts(1,3), '\color{blue}L_1');

plot3(rvInt(:,1), rvInt(:,2), rvInt(:,3), 'r', 'LineWidth', 1.2);
plot3(rvTgt(idx_f,1), rvTgt(idx_f,2), rvTgt(idx_f,3), '.r', 'MarkerSize', 14);
plot3(rv(:,1), rv(:,2), rv(:,3), 'w', 'LineWidth', 1.5);
axis equal; grid on; set(gca,'Clipping','off');
```
# Notes and Troubleshooting

\- If `tfMin` fails to converge: \* Adjust the costate guess (`lambda0_tf`) \* Use continuation on `Tmax` or final time \* Reduce solver tolerances or switch ODE integrator \- Verify `pumpkyn.cr3bp.tfMin` and `.tfMinProp` signatures match your toolbox version.

