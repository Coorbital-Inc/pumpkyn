function haloTransferDemo()

%% Characteristic mass (sum of primaries)
                          M = 5.9736E24 + 7.35E22;  

%% Initial orbit from interpolating seeds:
                        Np = 1;       %L2 Halo Orbit
                      tau0 = 1.085*pi; 
                        pm = +1;
[tau0,x0,muStar,tStar,lStar] = pumpkyn.cr3bp.getTulip(tau0,Np,pm,1e-8);
%L1 halo orbit
%                          x0 = [8.2338739411138306E-1	0 6.9135470535882995E-3	0	1.2711946779802655E-1	0];
%                        tau0 = 2.7433219717780992E+0;
%                      muStar = 0.012150585609624;
%                       lStar = 389703.264829278;
%                       tStar = 382981.289129055;
              [tauInt,rvInt] = pumpkyn.cr3bp.prop(tau0,x0,muStar);

%% Desired orbit to transfer to from interpolating seeds:
                      Np = 7; %Number of petals
                      pm = -1;     
                    tau0 = (5/6)*2*pi;
                  [~,x0] = pumpkyn.cr3bp.getTulip(tau0,Np,pm,1e-8);
          [tauTgt,rvTgt] = pumpkyn.cr3bp.prop(tau0,x0,muStar);


%% Pick a point on the halo orbit and compute a lambert transfer:
pumpkyn.cr3bp.showMoon(lStar,muStar);
plot3(rvInt(:,1),rvInt(:,2),rvInt(:,3),'b');
plot3(rvTgt(:,1),rvTgt(:,2),rvTgt(:,3),'y');
axis equal;

                    tof = 3.5*86400/tStar;
                 tauVec = linspace(0,tauTgt(end) - tauTgt(end)/21,21)';
                tauVec0 = linspace(0,tauInt(end) - tauInt(end)/9,9)';

                    
                  dVtot = []; 
                 
  for ti=1:numel(tauVec0)   
                        rv0 = interp1(tauInt,rvInt,tauVec0(ti),'spline');
                        plot3(rv0(1,1),rv0(1,2),rv0(1,3),'.g','markersize',20);
       v0tmp = [];
       dVtmp = [];

      for tt=1:numel(tauVec)               
                        rvf = interp1(tauTgt,rvTgt,tauVec(tt),'spline');
          [v0,vf,converged] = pumpkyn.cr3bp.minLambert(rv0,rvf,tof, ...
                                                       muStar,tStar,lStar,M,2);
        % Propagate & Plot the transfer trajectory:
        if converged
              dV1 = v0 - rv0(:,4:6);
              dV2 = rvf(:,4:6) - vf;
            dVtot = [dVtot; norm(dV1)+norm(dV2)];
            dVtmp = [dVtmp; dVtot(end)];
            v0tmp = [v0tmp; v0];
        end

      end
      %Show the min dV transfer:
               [~,idx] = min(dVtmp);
      [tauXfer,rvXfer] = pumpkyn.cr3bp.prop(tof,[rv0(:,1:3),v0tmp(idx,:)],muStar);
      plot3(rvXfer(:,1),rvXfer(:,2),rvXfer(:,3),'-','color',[1 1 1 0.4]);
      plot3(rvXfer(end,1),rvXfer(end,2),rvXfer(end,3),'.r','markersize',15);
      drawnow;

  end

  title(['\color{white}Minimum \DeltaV_{tot} = ',num2str(min(dVtot)*1e3*lStar/tStar),' [m/s]']);

end