function lowThrustLEO2GEO()
%% Purpose:
%
%  This routine will use the pumpkyn toolkit to solve for the low thrust
%  trajectory for a satellite to transfer from a LEO orbit to a GEO orbit
%  using the minimum time CR3BP formulation.
%

%% Constants:
  muStar = 0.012150585609624;          % Mass ratio
   lStar = 389703.264829278;
   tStar = 382981.289129055;
       M = 5.9736E24 + 7.35E22;        % Combined mass of primaries (kg)

%% Initial Orbit Conditions:
      a0 = 7000;         %semi-major axis (km)
      e0 = 0.0;          %eccentricity
      i0 = 28.5.*pi/180; %inclination (rad)
      O0 = 0;            %RAAN (rad)
      o0 = 0;            %ARGP (rad)
      M0 = mod(-220.*pi/180,2*pi); %Mean Anomaly (rad)
     nu0 = computeTrueAnomaly(e0,M0);

%% Final Orbit Conditions:
      aF = 42000;
      eF = 1e-3;
      iF = 1*pi/180;
      OF = 0;
      oF = 0;
      MF = 207.155823*pi/180;  %Optimized departure and arrival times
     nuF = computeTrueAnomaly(eF,MF);
 
%% Convert States from Keplerian -> PCI -> CR3BP:
  rvND0 = pumpkyn.cr3bp.fromOrb(0,[a0,e0,i0,o0,O0,nu0],M,muStar,tStar,lStar,1);
  rvNDF = pumpkyn.cr3bp.fromOrb(0,[aF,eF,iF,oF,OF,nuF],M,muStar,tStar,lStar,1);
  
%% Run Solver:
  m0 = 150;                                    % Initial mass (kg)
  g0 = 9.80665 * tStar^2 / (1000 * lStar);     % Gravity in ND units
aMax = 0.715e-01;                              % Avg Acceleration (m/s^2)
Tmax = m0*aMax;                                % Max thrust (N)
Tmax = (Tmax / m0) * tStar^2 / (lStar * 1000); % ND thrust accel
 Isp = 900 / tStar;                            % ND specific impulse
   c = Isp * g0;                               % ND exhaust velocity

%%  Solve Minimum-Time Transfer
% Use indirect optimal control to find the minimum-time trajectory.
    lambda0_tf = [ 1.68103229671784
                  -0.835819604294616
                   0.29209265591078
                   0.00225127896189283
                   0.00117398737057382
                  -0.000856414909752053
                   0.110719156572044
                   0.157957584644197];
         
sol_lambda0_tf = pumpkyn.cr3bp.tfMin(rvND0, rvNDF, lambda0_tf, Tmax, c, muStar); %,'lsqnonlin');

%% Propagate the Solution:
[tau, rv] = pumpkyn.cr3bp.tfMinProp(sol_lambda0_tf(8), ...
              [rvND0, 1, sol_lambda0_tf(1:7)'], Tmax, c, muStar);

%% Extract the Acceleration Due to Thrust and Hamiltonian of System:
rvDot = NaN(size(rv));
    H = NaN(size(tau));
   aT = NaN(numel(tau),3);
for tt=1:numel(tau)
    [rvDot(tt,:), H(tt),~,aT(tt,:)] = pumpkyn.cr3bp.tfMinEoM(tau(tt),rv(tt,:)',Tmax,c,muStar);
end

%Dimensionalize the acceleration due to thrust:
    aT_kms = aT.*lStar./(tStar.^2);
aT_Mag_kms = pumpkyn.util.vmag(aT_kms,2);
dV_tot_kms = trapz(tau.*tStar,aT_Mag_kms,1);


%% Show Pitch and Yaw Angles:
[~,theta,psi] = pumpkyn.cr3bp.toModEq(tau,rv,aT,M,muStar,tStar,lStar,1);
figure('color',[1 1 1]);
yyaxis left
plot(tau.*tStar./3600,theta.*180/pi,'k','linewidth',1.5);
grid on; xlabel('Time [Hr]'); ylabel('Pitch Angle [deg]');
set(gca,'YColor','k')    % left axis color
yyaxis right
plot(tau.*tStar./3600,psi.*180/pi,'r','linewidth',1.5);
grid on; xlabel('Time [Hr]'); ylabel('Yaw Angle [deg]');
set(gca,'YColor','r')    % left axis color

%% Show Pos/Vel costates:
figure('Color',[1 1 1]);
subplot(1,3,1)
plot(tau.*tStar./3600,rv(:,8:10),'-','linewidth',1.5);
grid on;
xlabel('Time [Hrs]'); ylabel('Position Costates');
legend('\lambda_x','\lambda_y','\lambda_z');
subplot(1,3,2)
plot(tau.*tStar./3600,rv(:,11:13),'-','linewidth',1.5);
grid on;
xlabel('Time [Hrs]'); ylabel('Velocity Costates');
legend('\lambda^{\prime}_x','\lambda^{\prime}_y','\lambda^{\prime}_z');
subplot(1,3,3)
plot(tau.*tStar./3600,rv(:,14),'-k','linewidth',1.5);
grid on;
xlabel('Time [Hrs]'); ylabel('Mass Costate');

%% Propellant Usage and ΔV Analysis
% Compute propellant consumption and cumulative ΔV.
          mTot = rv(:,7) .* m0;                      % Mass vs. time (kg)
         mProp = m0 - mTot;                          % Propellant used (kg)
         dVtot = (c.*lStar./tStar).*log(m0 ./ mTot); % ΔV (km/s)
         
%% Show 3D Trajectory in CR3BP:
hIn = figure('color', [0 0 0]);
pumpkyn.cr3bp.showEarth([],lStar, muStar, hIn);
plot3(rv(:,1),rv(:,2),rv(:,3),'w');
quiver3(rv(:,1),rv(:,2),rv(:,3),aT(:,1),aT(:,2),aT(:,3),'b');
plot3(rvND0(:,1),rvND0(:,2),rvND0(:,3),'.g','markersize',12);
plot3(rvNDF(:,1),rvNDF(:,2),rvNDF(:,3),'.r','markersize',12);
axis equal off;
set(gca, 'color', 'k','clipping','off');

%% Show Orbital Elements:
oev = pumpkyn.cr3bp.toOrb(tau,rv,M,muStar,tStar,lStar,1);
figure('color',[1 1 1]);
% Left axis for semi-major axis
yyaxis left
plot(tau .* tStar ./ 3600, oev(:,1) ./ 1000, 'k', 'linewidth', 2); hold on;
[aMax,idxMax] = max(oev(:,1));
text(tau(idxMax).*tStar./3600, (aMax/1000)*1.02, 'a', ...
     'Color','k','FontSize',12,'FontWeight','bold');
ylabel('Semi-major Axis [Mm]');
grid on
set(gca,'YColor','k')    % left axis color
% Right axis for eccentricity
yyaxis right
plot(tau .* tStar ./ 3600, oev(:,2), 'k', 'linewidth', 1.5); hold on;
[eMax,idxMax] = max(oev(:,2));
text(tau(idxMax).*tStar./3600, (eMax)*1.02, 'e', ...
     'Color','k','FontSize',12,'FontWeight','bold');
set(gca,'YColor','k')    % right axis color
ylabel('Eccentricity'); xlabel('Time [Hrs]'); xlim([0 17]);
title('Figure 1: Semi-major axis and ecc for a LEO-GEO Transfer');

figure('color',[1 1 1]);
% Left axis for semi-major axis
yyaxis left
plot(oev(:,1)./1000, oev(:,3).*180/pi, 'k', 'linewidth', 2); hold on;
[iMax,idxMax] = max(oev(:,3));
text(oev(idxMax,1)./1000, (iMax*180/pi)*1.02, 'i', ...
     'Color','k','FontSize',12,'FontWeight','bold');
ylabel('Inclination [deg]'); grid on; set(gca,'YColor','k')  
% Right axis for eccentricity
yyaxis right
plot(oev(:,1)./1000, oev(:,2), 'k', 'linewidth', 1.5); hold on;
[eMax,idxMax] = max(oev(:,2));
text(oev(idxMax,1)./1000, (eMax)*1.02, 'e', ...
     'Color','k','FontSize',12,'FontWeight','bold');
set(gca,'YColor','k');    % right axis color
ylabel('Eccentricity'); grid on; xlabel('Semi-Major Axis [Mm]');
title('Figure 2: Ecc and Inc vs semi-major axis for a LEO-GEO Transfer');

%% Show Delta-V:
figure('Color',[1 1 1]);
yyaxis left
[dVtotMax,idxMax] = max(dVtot);
plot(tau.*tStar./3600,dVtot,'k','linewidth',1.5);
grid on; ylabel('\DeltaV [km/s]');
set(gca,'YColor','k')    % left axis color
text(tau(idxMax,1).*tStar./3600, dVtotMax*1.01, '\DeltaV', ...
     'Color','k','FontSize',12,'FontWeight','bold');
yyaxis right
[aTtMax,idxMax] = max(aT_Mag_kms);
plot(tau.*tStar./3600,aT_Mag_kms.*1000*100,'k','linewidth',1.5);
grid on; ylabel('|a_T| [cm/s^2]'); xlabel('Time [Hrs]');
set(gca,'YColor','k')    % right axis color
aAvg = dVtot(end)/(tau(end)*tStar);
text(tau(idxMax,1).*tStar./3600, aTtMax*1.01.*1000.*100, '|a_T|', ...
     'Color','k','FontSize',12,'FontWeight','bold');
title(['Average Acceleration = ',num2str(aAvg*1000*100),' [cm/s^2]']);

end

function nu0 = computeTrueAnomaly(e,M)
   e2 = e.*e;
   e3 = e2.*e;
%Compute the guess for True anomaly:
  nu0 = M + (2.*e - e3./4).*sin(M) + ...
             (5.*e2.*sin(2.*M)./4) + ...
             (13.*e3.*sin(3.*M)./12);
%Mean anomaly can be computed from the eccentricity and true-anomaly
% See https://en.wikipedia.org/wiki/Mean_anomaly
 nu0 = fzero(@(nu)computeMeanAnomaly(e,nu) - M, nu0);
end

function M = computeMeanAnomaly(e,nu)
         M = atan2(sqrt(1-e.^2).*sin(nu),e + cos(nu)) - ...
             e.*sqrt(1 - e.^2).*sin(nu)./(1 + e.*cos(nu));
         M = mod(M,2*pi);
end