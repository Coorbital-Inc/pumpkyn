%% Tulip-Shaped Orbit Pseudo-Arclength Continuation
%
% This tutorial demonstrates how to extend a tulip-shaped orbit family 
% using pseudo-arclength continuation in the Earth–Moon CR3BP.
%
% The continuation process incrementally adjusts an existing periodic orbit 
% (the "seed") to find nearby family members with slightly different 
% energy or period. Repeating this procedure allows tracing an entire 
% branch of orbits until a desired target period is reached.
%
% This is a key technique for exploring resonant families and understanding 
% how orbit geometry evolves as a function of energy.

%% Specify the Tulip-Family:
                 Np = 7;          %Number of petals
                 pm = +1;         %Prograde or retrograde
% Retrieve an initial "seed" orbit for the specified family.
% Inputs:
%   tau0 – initial period (0 = use default catalog value)
%   Np   – number of petals
%   pm   – hemisphere selector
%
% Outputs:
%   tau0  – dimensionless orbital period of the seed orbit
%   mu    – Earth–Moon mass ratio
%   lStar – characteristic length scale [km]
[tau0,~,mu,~,lStar] = pumpkyn.cr3bp.getTulip(0,Np,pm);

%% Run Pseudo-Arclength Continuation Until a Sidereal Period is Reached
%
% The pseudo-arclength continuation method incrementally varies the period (tau0)
% and corrects the initial state (x0) to converge on a nearby periodic orbit.
% This allows the orbit family to be traced continuously, even through
% bifurcation or turning points where traditional continuation fails.

% Define the colormap used to color orbits as the period increases
 hsvVec = jet;
% Define the target maximum period (2π is a lunar sidereal month)
tau0Max = 2*pi;
% Generate a color mapping that corresponds to each intermediate period
 cVec = linspace(3.14,3.19,size(hsvVec,1))'; %linspace(tau0,tau0Max,size(hsvVec,1))';
% Initialize step size (pseudo-arclength parameter) and iteration counter
     ds = 0;
  count = 0;
 nullDF = zeros(7,1); 

while tau0 < tau0Max     
    if count == 0
        pumpkyn.cr3bp.showMoon(lStar,mu); hold on;
    end
    % Retrieve the initial conditions for the current tulip family member.
    % The period (tau0) is updated each loop iteration.
   [tau0,x0,mu] = pumpkyn.cr3bp.getTulip(tau0,Np,pm);
   % Perform pseudo-arclength continuation:
   %   - Perturbs the current orbit along the family tangent direction
   %   - Uses a correction loop to converge to a nearby periodic orbit
   %   - Returns the new state vector (x0) and corresponding updated period (tau0)
   %   - nullDF provides the tangent direction from the previous iteration
   [x0,tau0,~,nullDF] = pumpkyn.cr3bp.cont(x0(:),tau0,ds,nullDF,mu,1e-6);
                [~,x] = pumpkyn.cr3bp.prop([0 tau0],x0,mu);
                    C = pumpkyn.cr3bp.jacobi(x0(:)',mu);
   fprintf('Computing Period %f of %f, C = %f\n',tau0,tau0Max,C);
   plot3(x(:,1),x(:,2),x(:,3),'color',interp1(cVec,hsvVec,C,'linear')); 
   drawnow;
   % Increment the period and counter for the next continuation step
   tau0 = tau0+0.075;
   count = count+1;
end
colormap(hsvVec); clim([min(cVec), max(cVec)]); 
hC = colorbar; hC.Title.String ='\color{white}Jacobi Constant'; hC.Color = 'w';
% Interpretation:
%   - Each curve represents a unique member of the 7-petal tulip family.
%   - As the period increases, the orbit's geometry and energy change smoothly.
%   - The color gradient encodes the orbit's relative energy level.
%   - Pseudo-arclength continuation enables robust exploration of orbit families
%     even through turning points where simple parameter continuation fails.