%% Converting a Pumpkin Orbit from CR3BP to J2K Coordinates
%
% This tutorial demonstrates how to convert a *pumpkin-shaped* orbit
% defined in the **barycentric rotating frame** of the
% Earth–Moon Circular Restricted Three-Body Problem (CR3BP)
% into the **inertial J2000 (J2K)** frame.
%
% Such conversions allow us to visualize and analyze CR3BP orbits
% in an inertial reference frame used by standard mission design tools.
% This enables higher-fidelity modeling, comparison with ephemeris data,
% and low-energy transfer trajectory studies.

%% Establish physical and system constants
jd0 = juliandate(2020,1,1);     % Initial epoch (Julian date)
G   = 6.67384e-20;              % Universal gravitational constant [km^3/kg/s^2]
M   = 5.9736E24 + 7.35E22;      % Combined mass of Earth and Moon [kg]

% The CR3BP formulation models the Earth–Moon system as two massive primaries
% (Earth and Moon) moving in circular orbits about their barycenter.
% The spacecraft's mass is assumed negligible (the "restricted" case).

%% Define a 13-petal pumpkin-shaped orbit in the CR3BP rotating frame
tau0 = 2*pi;     % Dimensionless orbital period (1 rotation in normalized time)
Np   = 13;       % Number of petals in the pumpkin orbit family
pm   = +1;       % +1 = Northern family, -1 = Southern family

% Retrieve initial state and scaling parameters for this orbit family
[tau0, x0, muStar, tStar, lStar] = pumpkyn.cr3bp.getTulip(tau0, Np, pm, 1e-12);

%% Propagate the orbit in the CR3BP rotating frame
% This integrates the nonlinear CR3BP equations of motion for one full period.
[tau, x] = pumpkyn.cr3bp.prop(tau0, x0, muStar);

%% Convert from CR3BP rotating frame to inertial (J2K) frame
% The resulting inertial position and velocity vectors allow direct comparison
% with traditional mission analysis outputs (e.g., STK, GMAT, or SPICE).
jd = jd0 + tau .* tStar ./ 86400;    % Convert dimensionless time → Julian Date
rv = pumpkyn.cr3bp.toJ2K(jd, x, muStar, M, 2);   % Convert to J2000 frame

%% Visualize the trajectory in inertial space
% Plot the Moon as a 3D sphere and overlay the converted pumpkin orbit.
pumpkyn.util.moon3D();   % Draw Moon surface in J2000 coordinates
plot3(rv(:,1), rv(:,2), rv(:,3), '-', 'Color', [0.8500, 0.3250, 0.0980 1], ...
'LineWidth', 1); axis equal;
set(gca, 'Clipping', 'off','view', [35.0989 0]);
title('\color{orange}Inertial J2000 Frame');
