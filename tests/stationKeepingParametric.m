function stationKeepingParametric()

figure('Color',[1 1 1]);


                          Np = 2:20; %Number of petals
                        tau0 = 2*pi; 
                          pm = +1;
                          dt = (1:0.1:20)'.*86400; %seconds 

          tiledlayout(4,5,'TileSpacing','tight');

          for tp=1:numel(Np)
              nexttile; 
              [dV_LB,dV_UB] = computeSKParametric(Np(tp),tau0,pm,dt);
              plot(dt./3600,dV_LB,'b'); hold on;
              plot(dt./3600,dV_UB,'r');
              grid on;
              xlabel('\Deltat_b [Hrs]'); ylabel('\DeltaV_{tot} [m/s/yr]');
              title(['N_p = ',num2str(Np(tp))]);
              drawnow;
          end


end

function [dV_LB,dV_UB] = computeSKParametric(Np,tau0,pm,dt)
[tau0,x0,muStar,tStar,lStar] = pumpkyn.cr3bp.getTulip(tau0,Np,pm,1e-12);
                    dVMaxErr = sqrt(3)*(0.01/1000)*tStar/lStar; % cm/s -> ND
                        dTau = dt./tStar; %Correction interval Days -> ND
                           N = (365.25*86400/tStar)./dTau;
                       dV_LB = NaN(size(N));
                       dV_UB = NaN(size(N));
          for tt=1:numel(dTau)            
            [dV_LB(tt,1),dV_UB(tt,1)]  = pumpkyn.cr3bp.stationKeeping_deltaV(x0,tau0, ...
                                         dVMaxErr,dTau(tt),muStar,N(tt));
          end
          dV_LB = (dV_LB.*lStar./tStar).*1000; %m/s
          dV_UB = (dV_UB.*lStar./tStar).*1000; %m/s   
end