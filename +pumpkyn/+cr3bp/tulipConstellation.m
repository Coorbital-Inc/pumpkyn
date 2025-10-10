function [tau,r,v,mu,lStar,tStar] = tulipConstellation(Np,tau0,Nr,pm,dtau)
%% Purpose:
%
%  This routine will generate a constellation of satellites spaced equally
%  in time using the CR3BP equations of motion, it will then output
%  the position and velocities in Non-Dimensional rotating barycenter
%  frame of reference.
%
%% Inputs:
%
%  Np               Integer                     Number of Petals Used
%                                               for tulip-shaped orbit
%
%  tau0             double                      Dimensionless period
%                                               tau0 = q*2*pi/p;
%
%  Nr               Integer                     Number of revolutions to
%                                               propagate the CR3BP orbit
%                                               tauMax = Nr*tau0
%
%  pm               Integer                     +1 = Northern Hemisphere
%                                               -1 = Southern Hemisphere
%
%  dtau             [1 x Ns]                    Offset vector corresponding
%                                               to each satellite location
%                                               (dimensionless)
%
%% Outputs:
%
%  tau              [N x 1]                     Time in from
%                                               epoch (dimensionless)
%
%  r                [N x 3 x Ns]                Position of each satellite
%                                               (dimensionless)
%
%  v                [N x 3 x Ns]                Velocity of each satellite
%                                               (dimensionless)
%
%% Revision History:
%  Darin C. Koblick                                              07/30/2025
%  Copyright 2025 Coorbital, Inc.
%% ------------------------ Begin Code Sequence ---------------------------
if nargin == 0
                    %tStar = 382981.289129055;
                       Np = 5;
                        p = 5;
                        q = 4;
                       Nr = 2;
                       pm = -1;
                     tau0 = q*2*pi/p;
                       Ns = Np;
                     dtau = (linspace(0,Ns-1,Ns)).*tau0/Np;
       [tau,r,~,mu,lStar] = pumpkyn.cr3bp.tulipConstellation(Np,tau0,Nr,pm,dtau);
     
     %% Show the Moon Centered Inertial trajectory:
      %r = r - [1-mu,0,0];
     pumpkyn.cr3bp.showMoon(lStar,mu);
     plot3(r(:,1,1),r(:,2,1),r(:,3,1),'w'); 
     %Interpolate position to ensure equal time spacing:
     tau_i = linspace(0,tau(end),Nr*Np*90)';
         r = interp1(tau,r,tau_i,'spline');
     for ts=1:Ns  
         h{ts,1} = plot3(r(1,1,ts), ...
                         r(1,2,ts), ...
                         r(1,3,ts),'.','markersize',15);
     end
     axis equal;
     for tt=1:size(r,1)
        for ts=1:Ns
            set(h{ts,1},'XData',r(tt,1,ts), ...
                        'YData',r(tt,2,ts), ...
                        'ZData',r(tt,3,ts));
        end
        drawnow;
     end
    [tau,r,v,mu,lStar,tStar] = deal([]);
    return;
end

%% Lookup Orbit from Seed States:
[tau0,x0,mu,tStar,lStar] = pumpkyn.cr3bp.getTulip(tau0,Np,pm);
          
  
%% Propagate the Orbit for first satellite (only 1 revolution):   
          [tauND,rvND] = pumpkyn.cr3bp.prop(tau0,x0,mu);
                 tau_i = tauND + dtau;
                rvND_i = interp1(tauND,rvND,mod(tau_i,tau0),'spline');
                if size(tau_i,2) > 1
                    rvND_i = permute(rvND_i,[1 3 2]);
                end
                  
%% Replicate Solution:
                   tau = tauND;
                  rvND = rvND_i;
for ts=2:Nr
                   tau =  [tau(1:end-1); + tau(end)+tauND];
                   rvND = [rvND(1:end-1,:,:); rvND_i];
end

%% Outputs:
              r = rvND(:,1:3,:);
              v = rvND(:,4:6,:);
end



