function propPumpkin()
%% Retreive a tulip-shaped Orbit seed:
% Define the desired tulip family parameters:
                        tau0 = 2*pi;
                          Np = 20;
                          pm = +1;
% Retrieve the initial state vector and system constants.
% Inputs:
%   tau0 – Desired dimensionless orbital period
%   Np   – Number of petals in the family
%   pm   – Hemisphere flag
%   tol  – Continuation tolerance for orbit closure
%
% Outputs:
%   x0    – Dimensionless initial state vector [x, y, z, xdot, ydot, zdot]
%   mu    – Mass ratio of the Earth–Moon system
%   tStar – Characteristic time scaling [s]
%   lStar – Characteristic length scaling [km]
[tau0, x0, mu, tStar, lStar] = pumpkyn.cr3bp.getTulip(tau0,Np,pm,1e-12);
% Propagate one full period of the orbit using the CR3BP equations of motion.
                    [tau, x] = pumpkyn.cr3bp.prop(tau0,x0,mu);

%% Plot Orbit in the Rotating Frame
%
% This section visualizes the pumpkin-shaped orbit in the barycentric
% rotating frame, with the Moon shown for reference. The pumpkin orbit
% exhibits repeating loops (petals) around the Moon, forming a distinct
% 13-petal pattern when viewed over a full cycle.

% Plot the Moon in barycentric coordinates at its equilibrium location [1 - mu, 0, 0].
 pumpkyn.cr3bp.showMoon(lStar,mu);
% Plot the propagated orbit in orange to emphasize its geometry.
 plot3(x(:,1),x(:,2),x(:,3),'-', ...
      'Color', [[255, 92, 0]./256,0.2],'LineWidth',3);
 axis equal; set(gca,'clipping','off'); view(290,0);
 title(['\color{orange}',num2str(Np),'-Petal Pumpkin Orbit']);

%% Compute and Display Orbit Properties
%
% Compute and display key physical and dynamical properties of the orbit,
% including its stability index, energy (Jacobi constant), and altitudes.
% These values help characterize the orbit's behavior and mission relevance.
 data = pumpkyn.cr3bp.orbitProperties(x0,tau0,mu,lStar);
rMoon = 1737.1;
disp(['Stability index: ', num2str(data.StabilityIndex)]);
disp(['Jacobi Constant: ', num2str(data.Jacobi)]);
disp(['Perilune Alt: ',    num2str(data.Perilune.*lStar - rMoon),' [km]']);
disp(['Apolune Alt:  ',    num2str(data.Apolune.*lStar - rMoon),' [km]']);
disp(['Total Lunar Occultation:  ',    num2str(data.TotLunarOcc.*tStar./3600),' [Hrs]']);

end

