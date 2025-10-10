%% Station Keeping Estimation
% Compute the upper and lower bounds required to station keep in a pumpkin
% orbit.
%
% This tutorial demonstrates:
%   1. Generating a pumpkin orbit
%   2. Setting up station keeping constraints
%   3. Plotting and visualizing bounds
%
% Station-keeping analysis is essential for understanding how much ΔV is 
% required to correct small deviations from a periodic orbit due to model 
% inaccuracies, navigation errors, or dynamical instabilities.

%% Extract pumpkin orbit:
% Define the desired pumpkin orbit parameters:
                          Np = 13; 
                        tau0 = 2*pi; 
                          pm = +1;
% Retrieve the initial orbit state and system constants.
% Outputs:
%   x0    – Initial state vector [x, y, z, xdot, ydot, zdot]
%   muStar– Earth–Moon mass ratio (m2 / (m1 + m2))
%   tStar – Characteristic time scaling [s]
%   lStar – Characteristic length scaling [km]                          
[tau0,x0,muStar,tStar,lStar] = pumpkyn.cr3bp.getTulip(tau0,Np,pm,1e-12);

%% Set Up the Station-Keeping Constraints
%
% The following parameters define the correction schedule and allowable
% velocity dispersions used in the ΔV estimation.

% Define the maximum velocity error (dispersion) to be corrected per burn.
% A 2 cm/s dispersion is converted to dimensionless CR3BP units.
                dVMaxErr = sqrt(3)*(0.02/1000)*tStar/lStar; % 2 cm/s -> ND
% Define correction intervals from 1 to 60 days, converted to CR3BP units.
% Each value represents the time between station-keeping maneuvers.
                    dTau = (1:1:60)'.*86400/tStar;          % [ND]
% Compute the number of correction maneuvers per year based on dTau.
                       N = (365.25*86400/tStar)./dTau;      % # burns/year

%% Estimate Station-Keeping ΔV Bounds Using Linear Stability Theory
%
% The function pumpkyn.cr3bp.stationKeeping_deltaV() uses the eigenvalues of
% the monodromy matrix (linearized orbit dynamics) to estimate upper and 
% lower ΔV bounds. These correspond to:
%   - Lower Bound (dV_LB): Ideal correction, perfect control
%   - Upper Bound (dV_UB): Realistic correction, includes execution and 
%                          navigation errors (~factor of 3 margin)
dV_LB = NaN(size(N));
dV_UB = NaN(size(N));
for tt=1:numel(dTau)            
  [dV_LB(tt,1),dV_UB(tt,1)]  = pumpkyn.cr3bp.stationKeeping_deltaV(x0,tau0, ...
                                         dVMaxErr,dTau(tt),muStar,N(tt));
end
% Convert results from dimensionless units to physical units [m/s per year]
                       dV_LB = (dV_LB.*lStar./tStar).*1000; %m/s
                       dV_UB = (dV_UB.*lStar./tStar).*1000; %m/s    

%% Plot Station-Keeping Bounds
%
% This figure visualizes the total annual ΔV required to maintain a
% pumpkin-shaped orbit as a function of correction frequency.
% It shows the lower, upper, and average estimates, with uncertainty shading.
figure('Color',[1 1 1]);
dV_Avg = (dV_LB + dV_UB)./2;
    UB = dV_UB - dV_Avg;
    LB = dV_Avg - dV_LB;
% Plot shaded uncertainty region using a helper routine
pumpkyn.util.plotUnc(dTau.*tStar./86400,dV_Avg,UB,LB, ...
                     'color','k','faceAlpha',0.05,'EdgeColor','none');
hold on;
ax = gca;
ax.ColorOrder = [0 0 0; 1 0 0; 0 0 1];
plot(dTau.*tStar./86400,[dV_UB,dV_LB,dV_Avg],'LineWidth',2);
grid on;
xlabel('Maneuver Frequency [days]'); ylabel('\DeltaV_{tot} [m/s/yr]');
legend('','upper bound','lower bound','average value');
title(['Station-Keeping \DeltaV Estimation for N_p = ',num2str(Np)]);
% Interpretation:
%   - As correction frequency increases (shorter dTau), total ΔV decreases
%     because errors are corrected before they grow exponentially.
%   - The shaded region captures realistic uncertainty due to navigation 
%     and maneuver execution limits.