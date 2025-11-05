function [tau,y,Hf,dHdy] = tfMinProp(tf,y0,Tmax,c,muStar)
%% Purpose:
%
%  This routine will use event detection to stop and restart integration
%  when a switch point in tfMinEoM is detected.  This is the most accurate
%  and MATLAB-native way to detect a sign change in the switching function
%  S.
%
%% Inputs:
%
%  tf                       double                  Dimensionless Time
%                                                   corresponding to time
%                                                   of flight
%
%  y0                        [14 x 1]               Initial States
%                                                   Dimensionless:
%                                                   -------
%                                                   r [x,y,z]
%                                                   v [dr/dt]
%                                                   m [ND]
%                                                   lambda_r [x,y,z]
%                                                   lambda_v [x,y,z]
%                                                   lambda_m 
%
%  Tmax                     double                  Maximum Thrust
%                                                   Magnitude (N)
%
%  c                        double                  Exhaust Velocity
%                                                   (c = Isp*g0) of
%                                                   specific impulse
%
%  muStar                   double                  Mass Ratio of
%                                                   two primaries
%                                                   muStar = m2/(m1+m2)
%% Outputs:
%
%  tau                      [N  x 1]                Dimensionless Time
%
%  y                        [N x 14]                Pos/Vel/m/CoStates
%                                                   Dimensionless:
%                                                   -------
%                                                   r [x,y,z]
%                                                   v [dr/dt]
%                                                   m [ND]
%                                                   lambda_r [x,y,z]
%                                                   lambda_v [x,y,z]
%                                                   lambda_m 
%
% Hf                        double                  Hamiltonian at tf
%
%
% dHdy                      [1 x 14]                partial derivative
%                                                   of Hf wrt states
%
%% Revision History:
%  Darin C. Koblick                                          (c) 10/27/2025
%  Copyright Coorbital Inc.
%% ------------------------ Begin Code Sequence ---------------------------
if nargin == 0
   %Table Three: Minimum Time Solutions Computed:
   muStar = 1.21506683E-2;
    lStar = 3.84405000E5;   %characteristic length (km)
    tStar = 3.75676967E5;   %characteristic time (s)
       m0 = 1500;           %Initial Satellite mass (kg)
       g0 = 9.80665*tStar*tStar/(1000*lStar);       %m/s^2 -> ND
     Tmax = 10;                                     %Newtons
     Tmax = (Tmax/m0)*tStar^2/(lStar*1000);         %Max Thrust ND 
      Isp = 3000/tStar;    %Specific Impulse s -> ND
        c = Isp*g0;         %Exhuast Velocity (ND)

       %Boundary Conditions:
       r0 = [-0.019488511458668,-0.016033479812051,0];
       v0 = [8.918881923678198,-4.081793688818725,0];

       %Initial States:
       lambda_r0 = [15.616017, 32.875896, -0.094522];
       lambda_v0 = [-0.101606, 0.044791, -0.000150];
       lambda_m0 = 0.133266;
              tf = 8.6*86400/tStar;
              y0 = [r0(:); v0(:); m0/m0; lambda_r0(:); lambda_v0(:); lambda_m0];
 [tau,y,Hf,dHdy] = pumpkyn.cr3bp.tfMinProp(tf,y0,Tmax,c,muStar);

     figure('Color',[1 1 1]);
     plot3(y(:,1),y(:,2),y(:,3));
     axis equal;
     grid on;
     return;
end

     opts = odeset('RelTol',1e-10,'AbsTol',1e-12,'Events',@eventFunc);
        y = y0(:)';
      tau = 0;
    count = 0;
  while tau(end) < tf
        if count > 0
            %disp('switching');
        end
        [tau_t,y_t] = ode45(@pumpkyn.cr3bp.tfMinEoM, [tau(end) tf], y(end,:)',opts,Tmax,c,muStar);
                tau = [tau; tau_t(2:end)];      %#ok<AGROW>
                  y = [y; y_t(2:end,:)];        %#ok<AGROW>
              count = count+1;
  end

%% Re-compute the Final Hamiltonian and it's partial derivatives:
          yf = y(end,:)';
   if size(yf,1) < 210
         PHI = eye(14);
          yf = [yf; PHI(:)];
   end
 [~,Hf,dHdy] = pumpkyn.cr3bp.tfMinEoM(tau(end),yf,Tmax,c,muStar);
end


function [value, isterminal, direction] = eventFunc(~,y,~,c,~)
             m = y(7);
      lambda_v = y(11:13);
      lambda_m = y(14);
             S = -norm(lambda_v)*c/m - lambda_m;
         value = S;
    isterminal = 1;   % stop the integration
    direction  = 0;  % any direction
end

