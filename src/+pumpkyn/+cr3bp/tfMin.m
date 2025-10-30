function sol_lambda0_tf = tfMin(rv0,rvf,lambda0_tf,Tmax,c,muStar)
%% Purpose:
%
%   Solves a minimum-time transfer between two rotating-frame states in the
%   circular restricted three-body problem (CR3BP) using Pontryagin's
%   Minimum Principle (PMP) and a single-arc shooting method w an analytic
%   Jacobian derived from the propagated State Transition Matrix (STM).
%
%   The decision vector is x = [lambda0; tf], where lambda0 are the seven
%   initial costates (lambda_r[3], lambda_v[3], lambda_m) and tf is the
%   (dimensionless) final time. The solver enforces the terminal conditions:
%       r(tf) - r_f = 0,
%       v(tf) - v_f = 0,
%       lambda_m(tf) = 0,
%              H(tf) = 0,          (free-final-time transversality)
%   by solving R(x) = 0 with FSOLVE (trust-region-dogleg). The residual
%   Jacobian dR/dx is assembled analytically using the final STM and dH/dy.
%
%  ASSUMPTIONS / NOTES:
%
% • States are nondimensional CR3BP rotating-frame variables; mass is ND.
% • Control law uses min-time PMP (throttle & direction from switching fun).
% • For bang–bang solutions w switching, improve accuracy via event detect
%   and restarting integration @ switches.
% • Good initial guesses matter; thrust continuation (reduce Tmax gradually)
%   typically improves robustness.
%
%% Inputs:
%
%  rv0                      [1 x 6]                 Initial Pos/Vel @
%                                                   t = 0, dimensionless
%
%  rvf                      [1 x 6]                 Final Pos/Vel @
%                                                   t = tf, dimensionless
%
%  lambda0_tf               [8 x 1]                 Init guess of costates
%                                                   and final time of
%                                                   flight (dimensionless)
%                                                   lambda_r [x,y,z]
%                                                   lambda_v [x,y,z]
%                                                   lambda_m [double]
%                                                   tf       [double]
%
%
%  tf                       double                  Dimensionless Time
%                                                   corresponding to time
%                                                   of flight
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
%
%
%% Outputs:
%
% sol_lambda0_tf               [8 x 1]              solved costates
%                                                   and final time of
%                                                   flight (dimensionless)
%                                                   lambda_r [x,y,z]
%                                                   lambda_v [x,y,z]
%                                                   lambda_m [double]
%                                                   tf       [double]
%
%                                                   corresponds to:
%                                                   and rv(tf) - rvf = 0
%                                                       lambda_m(tf) = 0
%                                                              H(tf) = 0
%
%% Revision History:
%  Darin C. Koblick                                          (c) 10/29/2025
%  Copyright Coorbital Inc.
%% ------------------------ Begin Code Sequence ---------------------------

if nargin == 0
   useParallel = true;
   %Table Three: Minimum Time Solutions Computed:
   muStar = 1.21506683E-2;
    lStar = 3.84405000E5;                           %characteristic length (km)
    tStar = 3.75676967E5;                           %characteristic time (s)
       m0 = 1500;
       g0 = 9.80665 * tStar^2/(1000*lStar);         %m/s^2 -> ND
     Tmax = 0.8;                                    %Newtons
     Tmax = (Tmax/m0)*tStar^2/(lStar*1000);         %Max Thrust ND 
      Isp = 3000/tStar;    %Specific Impulse s -> ND
        c = Isp*g0;         %Exhuast Velocity (ND)
       %Boundary Conditions:
       rv0 = [-0.019488511458668,-0.016033479812051,0, ...
               8.918881923678198,-4.081793688818725,0];

       rvf = [0.823385182067467,0,-0.022277556273235, ...
              0,                0.134184170262437,0];

    %Initial States for Tmax = 0.8N:
       lambda0_tf = [ -78.2696859634265
         -380.485644187065
          5.76026987520663
          1.05579750843103
         -0.78118612362464
      -0.00743727273612572
          20.2103357875488
          15.3190846824972];

 sol_lambda0_tf = pumpkyn.cr3bp.tfMin(rv0,rvf, ...
                   lambda0_tf,Tmax,c,muStar,useParallel);

   %Take the solution and re-compute the trajectory: 
       [tau,rv] = pumpkyn.cr3bp.tfMinProp(sol_lambda0_tf(8), ...
                  [rv0,1,sol_lambda0_tf(1:7)'],Tmax,c,muStar);

     %Show the 3D trajectory:
      hIn = figure('color',[0 0 0]);
     pumpkyn.cr3bp.showEarth(lStar,muStar,hIn); hold on;
     pumpkyn.cr3bp.showMoon(lStar,muStar,hIn);
     plot3(rv(:,1),rv(:,2),rv(:,3),'w'); hold on;
     plot3(rv0(:,1),rv0(:,2),rv0(:,3),'.g','markersize',15);
     plot3(rvf(:,1),rvf(:,2),rvf(:,3),'.r','markersize',15);
     axis equal;
     grid on; set(gca,'Clipping','off');

     %Show the mass and delta-V:
     figure('Color',[1 1 1]);
     subplot(1,2,1);
     m = rv(:,7);
     plot(tau.*tStar./86400,m.*m0);
     grid on;
     xlabel('Time [Days]'); ylabel('Mass [kg]');
     subplot(1,2,2);
     dV = c.*log(1./m);
     plot(tau.*tStar./86400,dV.*lStar./tStar);
     grid on;
     xlabel('Time [Days]'); ylabel('\DeltaV [km/s]');

     return;
end

%Minimum Time Problem:
% Find (lambda_i, tf) such that [x(tf),lambda(tf)]' satisfies
%   r(tf) - rf = 0
%   v(tf) - vf = 0
% lambda_m(tf) = 0
%       Ht(tf) = 0
%
% Use a multi-variable shooting method to solve this problem, we have 8 
% equations and 8 unknowns.
%                                   
                  m0 = 1;
                opts = optimoptions('fsolve','Display','iter', ...
                                    'FunctionTolerance',sqrt(1E-12), ...
                                    'StepTolerance',1e-12, ...
                                    'FiniteDifferenceStepSize',1e-10, ...
                                    'MaxFunctionEvaluations',1e4, ...
                                    'SpecifyObjectiveGradient',true, ...
                                    'FiniteDifferenceType','central', ...
                                     'ScaleProblem', 'jacobian', ...
                                    'Algorithm','trust-region-dogleg');
%Use fsolve to determine the 7 costates and final time of flight:
[sol_lambda0_tf] = fsolve(@(x)shootingResidual(x,rv0(:),m0,rvf(:),Tmax,c,muStar), ...
                   lambda0_tf(:),opts);
end


function [R, J] = shootingResidual(x, rv0, m0, rvf, Tmax, c, muStar)
%----------------------------------------------------------------------
%  SHOOTINGRESIDUAL  Analytic residual and Jacobian for indirect CR3BP
%
%  Inputs:
%         x = [lambda0(7x1); tf]       initial costates + final time
%
%       rv0 = [1 x 6]                  initial pos/vel states
%
%        m0 = double                   initial mass
%
%       rvf = [1 x 6]                  final pos/vel states
%
%       Tmax = double                  Max acceleration
%
%          c = double                  Exhuast Velocity
%
%     muStar = double                  Mass ratio of primary bodies
%                                      muStar = m2/(m1+m2)
%
%  Outputs:
%     R  = residual vector for fsolve
%     J  = analytic Jacobian  dR/dx  (optional)
%----------------------------------------------------------------------
lambda0 = x(1:7);                          %Initial costates
     tf = x(8);                            %Final Time
% --- initial augmented state + STM ---
  y0_aug = [rv0(:); m0; lambda0(:)];

  if nargout > 1
    PHI0 = eye(14);
  y0_aug = [y0_aug; PHI0(:)];
  end
% --- integrate augmented EOM + STM ---
   [~,Y_aug] = pumpkyn.cr3bp.tfMinProp(tf,y0_aug,Tmax,c,muStar);
% --- extract final state and STM ---
        yf   = Y_aug(end,1:14).';
[Ff,Hf,dHdy] = pumpkyn.cr3bp.tfMinEoM(tf, yf, Tmax, c, muStar);

% --- build residual vector --------------------------------------------
R = [ yf(1:6) - rvf(:);                                  % Δr(tf) Δv(tf)
      yf(14);                                            %  λ_m(tf) = 0
      Hf ];                                              %  H(tf) = 0

    if nargout > 1
        %------------------------------------------------------------------
        %  Analytic Jacobian via STM
        %------------------------------------------------------------------
              PHI_f = reshape(Y_aug(end,15:210),14,14);
        % Derivatives wrt initial costates λ0 (columns 8:14)
        dxf_dlam0   = PHI_f(1:6, 8:14);                 % ∂(r,v)/∂λ0
        dlmm_dlam0  = PHI_f(14,  8:14);                 % ∂λ_m/∂λ0
        dH_dlam0    = dHdy * PHI_f(:, 8:14);            % ∂H/∂λ0   [1 x 7]
    
        % Derivatives wrt final time tf (chain rule: dy/dtf = F(yf))
        dxdtf   = Ff(1:6);                            % ∂(r,v)/∂tf
        dlmmdtf = Ff(14);                             % ∂λ_m/∂tf
        dHdtf   = dHdy*Ff;                            % ∂H/∂tf   [1 x 1]
    
        % Assemble full Jacobian dR/d[λ0, tf] [8 x 8]
        J = [ dxf_dlam0 , dxdtf;
              dlmm_dlam0, dlmmdtf;
              dH_dlam0  , dHdtf ];
    end

end