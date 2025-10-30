%% Low Thrust Minimum-Time Transfer: GEO -> Tulip-Shaped Orbit
% This live script demonstrates a minimum-time low-thrust transfer from
% a GEO-like equatorial orbit (Earth-centered) into a tulip-shaped
% orbit in the Earth–Moon CR3BP using the Pumpkyn toolbox.
%
% Sections:
%  1) Constants & nondimensional scaling
%  2) Define initial GEO state (J2000 / PCI -> CR3BP ND)
%  3) Fetch tulip seed and propagate both orbits
%  4) Select matching departure/arrival points
%  5) Setup low-thrust vehicle & solver params
%  6) Solve minimum-time problem (tfMin)
%  7) Propagate resulting low-thrust trajectory
%  8) Compute propellant & ΔV, and visualize
%
% Notes:
% - All CR3BP functions are called from `pumpkyn.cr3bp.*` API.
% - Variables with suffix `ND` are non-dimensional (CR3BP scaling).
% - Replace initial guesses for costates (lambda0_tf) with your own if needed.

%% 1) Characteristic quantities and CR3BP constants
% Characteristic mass (sum of primaries — Earth + Moon) and gravitational
% constant used for dimensional -> nondimensional conversions.
M = 5.9736E24 + 7.35E22;            % total primaries mass (kg)
G = 6.67384e-20;                    % gravitational constant (km^3 / kg / s^2) in chosen units
muStar = 0.012150585609624;         % mass ratio mu = m2/(m1+m2) (Earth-Moon)
lStar = 389703.264829278;           % characteristic length (km) used in Pumpkyn
tStar = 382981.289129055;           % characteristic time (s)

%% 2) Initial Low Earth Orbit (GEO) — dimensional and conversion to ND
% Define GEO radius (approx geostationary altitude) and convert to CR3BP ND.
muEarth = G*(1-muStar)*M;           % Earth's gravitational parameter used for GEO circular speed
rEarth = 6378;                      % Earth radius (km)
a0 = rEarth + 35786;                % GEO semi-major axis (km) — geostationary altitude above Earth

% Initial inertial state in an Earth-centered J2K / PCI frame:
% Position = [a0, 0, 0] km (equatorial)
% Velocity = [0, sqrt(muEarth/a0), 0] km/s (circular velocity)
% pumpkyn.util.Rx(0) is used for an identity rotation (placeholder for attitude)
rv0 = [a0, 0, 0, (pumpkyn.util.Rx(0) * [0, sqrt(muEarth/a0), 0]')'];

% Keplerian orbital period (seconds)
P0 = 2*pi*sqrt((a0^3)/muEarth);

% Convert the J2000/PCI GEO state to CR3BP nondimensional rotating barycenter frame
% Format: pumpkyn.cr3bp.fromPCI(jd0_or_flag, rv, mu, tStar, lStar, primary_flag)
% Here we call with "0" as a convenience (pumpkyn's signature accepts 0 for epoch if implemented).
rv0ND = pumpkyn.cr3bp.fromPCI(0, rv0, muStar, tStar, lStar, 1);

%% 3) Target Tulip orbit seed (lookup + natural-parameter cleanup)
% Choose tulip family parameters:
Np = 7;             % number of petals for tulip family
pm = -1;            % southern hemisphere (-1) or northern (+1)
tau0 = (5/6) * 2*pi; % dimensionless period guess for tulip (example)

% Get a seed for the tulip orbit (returns dimensionless state x0 for the requested tau0)
% pumpkyn.cr3bp.getTulip(tau0, Np, pm, tol) returns [~, x0] (or similar signature)
[~, x0] = pumpkyn.cr3bp.getTulip(tau0, Np, pm, 1e-12);

%% 4) Propagate both orbits to find sensible departure (GEO) and arrival (Tulip) times
% Propagate the tulip seed for one period (ND units)
[tauTgt, rvTgt] = pumpkyn.cr3bp.prop(tau0, x0, muStar);

% Propagate the nondimensionalized GEO orbit for one orbital period (converted to ND time)
[tauInt, rvInt] = pumpkyn.cr3bp.prop(P0 / tStar, rv0ND, muStar);

% Choose departure/arrival points by a simple heuristic: minimize radial distance or z (example)
% idx_f: index of tulip orbit point chosen as arrival (min of 5th column used as heuristic)
[~, idx_f] = min(rvTgt(:,5));  % adjust selection rule as required
[~, idx_0] = min(rvInt(:,1));  % departure selection on GEO propagation (example)

% Extract the chosen initial & final states (dimensionless)
rv0 = rvInt(idx_0, :);   % chosen GEO ND state at departure
rvf = rvTgt(idx_f, :);   % chosen tulip ND state at arrival

%% 5) Low-thrust vehicle properties and nondimensional conversions
m0 = 15;                                % initial spacecraft mass (kg)
g0 = 9.80665 * tStar^2 / (1000 * lStar);% convert g0 to ND units used by pumpkyn (ND acceleration units)
Tmax = 0.019;                           % Max thrust in Newtons (example)
% Convert Tmax (N) to nondimensional acceleration: (T/mass) scaled by tStar^2 / (lStar*1000)
Tmax = (Tmax / m0) * tStar^2 / (lStar * 1000);
Isp = 2100 / tStar;                     % Specific impulse converted to nondimensional seconds (Isp / tStar)
c = Isp * g0;                           % exhaust velocity in ND units (c = Isp * g0)

% NOTES:
% - pumpkyn.tfMin likely expects Tmax (ND acceleration), c (ND exhaust velocity),
%   and muStar. Confirm signature in your local pumpkyn implementation.

%% 6) Initial guess for costates / Lagrange multipliers (lambda0_tf)
% This initial guess is solver-sensitive. The vector contains costates for state
% variables + the unknown final time. Replace with your own if you have a better guess.
lambda0_tf = [ ...
   -5.0178526640353;
    3.17730030211899;
   -0.539375408679051;
    0.104084843259225;
   -0.200258747716391;
    0.0961706438135756;
    3.77874556494284;
    5.69600922675905
];

%% 7) Solve the minimum-time low-thrust transfer
% tfMin signature (example): sol = pumpkyn.cr3bp.tfMin(rv0, rvf, lambda0, Tmax, c, mu)
% - rv0: initial ND state (row, 6 elements)
% - rvf: final ND state (row, 6 elements)
% - lambda0_tf: initial guess for costates (vector)
% - Tmax: nondimensional max acceleration
% - c: nondimensional exhaust velocity
% - muStar: CR3BP mass ratio
%
% The solver returns a solution vector containing the optimized costates and final time.
sol_lambda0_tf = pumpkyn.cr3bp.tfMin(rv0, rvf, lambda0_tf, Tmax, c, muStar);

% Typical return: sol_lambda0_tf = [lambda1..lambda6, massAdj?, t_final] or similar.
% Check pumpkyn.tfMin docs to map indices. Here we assume:
% - sol_lambda0_tf(8) is the optimized final time (dimensionless)
% - sol_lambda0_tf(1:7) are the optimized costates/multipliers used for propagation

%% 8) Propagate the optimal low-thrust trajectory (tfMinProp)
% tfMinProp signature (example): [tau, rv] = pumpkyn.cr3bp.tfMinProp(t_final, initial_augmented_state, Tmax, c, mu)
% initial augmented state: [rv0, mass_fraction, lambda1..lambda6]'
% Here we propagate to reconstruct the transfer trajectory for plotting/analysis.
[tau, rv] = pumpkyn.cr3bp.tfMinProp(sol_lambda0_tf(8), ...
                  [rv0, 1, sol_lambda0_tf(1:7)'], Tmax, c, muStar);

% rv output expected columns (dimensionless): [x y z xdot ydot zdot massFraction ...]
% tau: nondimensional time vector from 0 to t_final

%% 9) Compute propellant used, ΔV profile, and plot results (dimensionalized)
figure('color', [0 0 0]);

% Total mass at each time (kg)
mTot = rv(:,7) .* m0;
mProp = m0 - mTot;

% Total delta-V (km/s) computed from ideal rocket equation:
% dV (ND) = c * ln(m0 / m(t)). Convert to km/s using lStar and tStar scaling.
dVtot = c .* log(m0 ./ mTot) .* lStar ./ tStar;

% Plot required prop mass vs transfer time (days)
subplot(1,2,1);
plot(tau .* tStar ./ 86400, mProp, 'w', 'linewidth', 2);
xlabel('Transfer Time [Days]');
ylabel('Required Prop Mass [kg]');
grid on;
set(gca, 'color', 'k', 'xcolor', 'w', 'ycolor', 'w');

% Plot ΔV profile vs transfer time (km/s)
subplot(1,2,2);
plot(tau .* tStar ./ 86400, dVtot, 'w', 'linewidth', 2);
xlabel('Transfer Time [Days]');
ylabel('\DeltaV_{tot} [km/s]');
grid on;
set(gca, 'color', 'k', 'xcolor', 'w', 'ycolor', 'w');

%% 10) 3D visualization of Earth, Moon, reference orbits, and the transfer
hIn = figure('color', [0 0 0]);
% showEarth/showMoon place the primary bodies scaled to lStar with correct positions
pumpkyn.cr3bp.showEarth(lStar, muStar, hIn);  % draws Earth (dimensionless scaled)
pumpkyn.cr3bp.showMoon(lStar, muStar, hIn);   % draws Moon at [1-mu, 0, 0]

set(gca, 'color', 'k');
hold on;

% Plot target tulip orbit (green) and the propagated GEO orbit (red)
plot3(rvTgt(:,1), rvTgt(:,2), rvTgt(:,3), 'g', 'linewidth', 1.2);
plot3(rvInt(:,1), rvInt(:,2), rvInt(:,3), 'r', 'linewidth', 1.2);

% Show selected arrival/departure nodes with markers
plot3(rvTgt(idx_f,1), rvTgt(idx_f,2), rvTgt(idx_f,3), '.r', 'markersize', 14);
plot3(rvTgt(idx_f+1,1), rvTgt(idx_f+1,2), rvTgt(idx_f+1,3), '.b', 'markersize', 12);
plot3(rvInt(idx_0,1), rvInt(idx_0,2), rvInt(idx_0,3), '.g', 'markersize', 14);
plot3(rvInt(idx_0+5,1), rvInt(idx_0+5,2), rvInt(idx_0+5,3), '.b', 'markersize', 12);

% Plot the propagated low-thrust transfer in white
plot3(rv(:,1), rv(:,2), rv(:,3), 'w', 'linewidth', 1.5);

grid on;
axis equal;
set(gca, 'clipping', 'off');

%% End of script
% Suggestions & troubleshooting:
% - If tfMin fails to converge, try:
%     * improving the initial lambda guess (lambda0_tf),
%     * using continuation on Tmax or final time,
%     * relaxing tolerances or switching integrators.
% - Verify pumpkyn.cr3bp.tfMin and pumpkyn.cr3bp.tfMinProp input/output signatures
%   if your local toolbox version differs from the assumptions here.