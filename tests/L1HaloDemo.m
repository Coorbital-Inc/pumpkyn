function L1HaloDemo()
%% L1 Halo Orbit Demo
% This script demonstrates how to generate and visualize a family of
% Earth–Moon L1 halo orbits using continuation methods implemented
% in the Pumpkyn CR3BP toolbox.
%
% The orbit is colored by Jacobi constant, and continuation proceeds until
% the orbit's perigee intersects the lunar surface.
%
% Author: Darin C. Koblick
% Coorbital, Inc. (c) 2025

%% Setup CR3BP Parameters
% Mass ratio of the Earth–Moon system:
           mu = 1.215058560962404E-2;

% Compute Lagrange points for reference:
         Lpts = pumpkyn.cr3bp.lagrangePts(mu);

% Characteristic scaling quantities:         
        lStar = 389703;        % Characteristic length (km)
        rMoon = 1737.1/lStar;  % Dimensionless lunar radius

%% Define Initial L1 Halo Seed
% Initial guess for a L1 halo orbit (dimensionless state):        
           x0 = [8.2339081983651485E-1 0 9.8941366235910004E-4 0  1.2634272983881797E-1 0];

%+Z = Northern Apolune, -Z = Southern Apolune
        x0(3) = -abs(x0(3));

% Dimensionless orbital period (tau0) of initial orbit:           
         tau0 = 2.7430007981241529E+0;

% Convergence tolerance and continuation parameters:
          tol = 1e-8;
       nullDF = zeros(7,1);
           ds = 0;

%% Refine the Initial Orbit
% Perform a single correction to ensure periodic closure:           
                  [x00,tau00] = pumpkyn.cr3bp.cont(x0,tau0,ds,nullDF,mu,tol);

% Extract useful orbit properties (Jacobi, trajectory, perilune, etc.)
                         data = pumpkyn.cr3bp.orbitProperties(x00(:)',tau00,mu,lStar);

%% Visualization Setup
% Prepare continuation step size parameters:   
   ds = -1e-3;  % Step direction
dsMax = 1e-1;  % Max step magnitude

% Create color map keyed by Jacobi constant:
cColor = jet();
jColor = linspace(2.95,3.18,size(cColor,1))';

% Create 3D figure with dark theme:
figure('Color',[0 0 0]);
plot3(data.x(:,1),data.x(:,2),data.x(:,3), ...
      '-','color',interp1(jColor,cColor,data.Jacobi,'nearest')); hold on;
set(gca,'color','k','xcolor','k','ycolor','k','Clipping','off');
colormap(cColor); % Apply the custom colormap
clim([jColor(1), jColor(end)]); % Set the colormap limits to match your data range
c = colorbar; % Get the colorbar handle
c.Label.String = '\color{white}Jacobi Constant'; % Set a label for the colorbar
c.Color = 'w';
       [x,y,z] = sphere(30);
       surf(x.*rMoon + (1 - mu), ...
            y.*rMoon, ...
            z.*rMoon,'FaceColor','w','edgeColor','none');

% Plot L2 points as white dots:       
plot3(Lpts(1,1),Lpts(1,2),Lpts(1,3),'.w','markersize',15);       
axis equal off;

%% Mirror Orbit about x-y plane
% flip the initial state from +z to -z and re-run continuation
      x0 = x00(:);
    tau0 = tau00;
  nullDF = zeros(7,1);

    %% Orbit Continuation Loop
    % Continue orbit family until perilune altitude becomes negative:
    while (data.Perilune - rMoon) > 0
                         xOld = x0;
                       tauOld = tau0;
        [x0,tau0,conv,nullDF] = pumpkyn.cr3bp.cont(x0,tau0,ds,nullDF,mu,tol);
                            J = pumpkyn.cr3bp.jacobi(x0(:)',mu);
        if conv &&  tau0 > 0.1 && J > 2.9 && J < 4
            ds = sign(ds)*min(abs(ds)*1.05,dsMax);
          data = pumpkyn.cr3bp.orbitProperties(x0(:)',tau0,mu,lStar);
          plot3(data.x(:,1),data.x(:,2),data.x(:,3), ...
                '-','color',interp1(jColor,cColor,J,'nearest'));
          drawnow;
        else
              ds = ds*0.5;
              x0 = xOld;
            tau0 = tauOld;
        end
        disp(J);
    end


%% Done!
% The continuation terminates when the orbit intersects the lunar surface.
% You now have a visual sweep of the L1 halo family with Jacobi gradient.
end