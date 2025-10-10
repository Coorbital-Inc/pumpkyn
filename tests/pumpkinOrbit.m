%% 
function pumpkinOrbit()
%% Purpose:
%
%  This routine will demonstrate how to show a pumpkin type orbit found
%  from performing continuation on the tulip-shaped orbit family.
%
%% Revision History:
%  Darin C. Koblick                                         (c) 09/30/2025
%% -------------------- Begin Code Sequence -------------------------------


recordMovie = false;

%% Get the tulip orbit from interpolating seeds:
                      Np = 14; %Number of petals
                    tau0 = 2*pi; 
                      pm = +1;
[tau0,x0,mu,tStar,lStar] = pumpkyn.cr3bp.getTulip(tau0,Np,pm);

%% Propagate the orbit for a full revolution:
            [tau,rv] = pumpkyn.cr3bp.prop(tau0,x0,mu);
               tau_i = linspace(0,tau0,2*360*Np)';
                  rv = interp1(tau,rv,tau_i,'spline');
                 tau = tau_i;

%% Place satellites in coordinates representing a mouth:
             rSmile = [0.9517,0.0136075,-0.00847807
                       0.950914,0.00267867,-0.0139984
                       0.951972,-0.00229963,-0.0166305
                       0.954445,-0.0134423,-0.016462
                       0.956186,-0.0131344,-0.0198675 %tooth
                       0.957977,-0.0167139,-0.0198905
                       0.956703,-0.0178243,-0.016462
                       0.96679,-0.032726,-0.00664261
                       0.963833,-0.0259917,-0.0175635
                       0.959973,-0.0151283,-0.0237488
                       0.95584,-0.00132106,-0.0237488
                       0.955051,0.0014451,-0.0222651
                       0.9517,0.0136075,-0.00847807];

             rEyeL = [0.951674,0.0114861,0.0111427
                      0.950288,0.00544934,0.0112024
                      0.963998,0.00553285,0.0313263
                      0.951674,0.0114861,0.0111427];

             rEyeR = [0.956398,-0.0210575,0.0114017
                      0.960432,-0.0262091,0.0110707
                      0.972171,-0.0128244,0.0341085
                      0.956398,-0.0210575,0.0114017]; 

%% Identify satellites which represent these orbits:
    [~,idx] = min(pumpkyn.util.vmag(rv(:, 1:3) - permute(rSmile,[3 2 1]),2));
      [idx] = ind2sub([size(rv,1),1,size(rSmile,1)],idx(:));
    rvSmile = rv(idx,:);

    [~,idx] = min(pumpkyn.util.vmag(rv(:, 1:3) - permute(rEyeL,[3 2 1]),2));
      [idx] = ind2sub([size(rv,1),1,size(rEyeL,1)],idx(:));
     rvEyeL = rv(idx,:);

    [~,idx] = min(pumpkyn.util.vmag(rv(:, 1:3) - permute(rEyeR,[3 2 1]),2));
      [idx] = ind2sub([size(rv,1),1,size(rEyeR,1)],idx(:));
     rvEyeR = rv(idx,:);

%% Backward Propagate states:
       dTau = -tau0/Np;
   rvSmile0 = propagateAllInitStates(rvSmile,dTau,mu);
    rvEyeL0 = propagateAllInitStates(rvEyeL,dTau,mu);
    rvEyeR0 = propagateAllInitStates(rvEyeR,dTau,mu);

%% Forward Propagate All States:
      tauPts = (0:Np)'.*tau0/Np; 
      tauVec = unique([tauPts; pumpkyn.util.chebyspace(0,tau0/Np,90)']);
     rvSmile = propagateAllStates(rvSmile0,tauVec,mu);
      rvEyeL = propagateAllStates(rvEyeL0,tauVec,mu);
      rvEyeR = propagateAllStates(rvEyeR0,tauVec,mu);

%% Show Propagated states:
            pumpkyn.cr3bp.showMoon(lStar,mu);
            plot3(rv(:,1),rv(:,2),rv(:,3),'-', ...
            'Color', [[255, 92, 0]./256,0.2],'LineWidth',3);
            axis equal;
            set(gca,'clipping','off')
            view(290,0);

            hSO = plot3(squeeze(rvSmile(1,1,:)),squeeze(rvSmile(1,2,:)),squeeze(rvSmile(1,3,:)), ...
                 '.','color','w','markersize',15);
            hLO = plot3(squeeze(rvEyeL(1,1,:)),squeeze(rvEyeL(1,2,:)),squeeze(rvEyeL(1,3,:)), ...
                 '.','color','w','markersize',15);
            hRO = plot3(squeeze(rvEyeR(1,1,:)),squeeze(rvEyeR(1,2,:)),squeeze(rvEyeR(1,3,:)), ...
                 '.','color','w','markersize',15);

             hSl = plot3(squeeze(rvSmile(1,1,:)), ...
                            squeeze(rvSmile(1,2,:)), ...
                            squeeze(rvSmile(1,3,:)), ...
                 '-','color',[1 1 1 0.5]);
                hLl = plot3(squeeze(rvEyeL(1,1,:)), ...
                            squeeze(rvEyeL(1,2,:)), ...
                            squeeze(rvEyeL(1,3,:)), ...
                 '-','color',[1 1 1 0.5]);
                hRl = plot3(squeeze(rvEyeR(1,1,:)), ...
                            squeeze(rvEyeR(1,2,:)), ...
                            squeeze(rvEyeR(1,3,:)), ...
                 '-','color',[1 1 1 0.5]);
            count = 0;

            if recordMovie
                v = VideoWriter('ex_pumpkin_orbit.avi','Motion JPEG AVI');
                v.Quality = 100;
                open(v);
            end

            for tt=1:size(rvSmile,1)
                count = count+1;
                set(hSO,'XData',squeeze(rvSmile(tt,1,:)), ...
                        'YData',squeeze(rvSmile(tt,2,:)), ...
                        'ZData',squeeze(rvSmile(tt,3,:)));
                set(hLO,'XData',squeeze(rvEyeL(tt,1,:)), ...
                        'YData',squeeze(rvEyeL(tt,2,:)), ...
                        'ZData',squeeze(rvEyeL(tt,3,:)));
                set(hRO,'XData',squeeze(rvEyeR(tt,1,:)), ...
                        'YData',squeeze(rvEyeR(tt,2,:)), ...
                        'ZData',squeeze(rvEyeR(tt,3,:)));

                if mod(tauVec(tt),tau0/Np) == 0
                set(hSl,'XData',squeeze(rvSmile(tt,1,:)), ...
                        'YData',squeeze(rvSmile(tt,2,:)), ...
                        'ZData',squeeze(rvSmile(tt,3,:)));
                set(hLl,'XData',squeeze(rvEyeL(tt,1,:)), ...
                        'YData',squeeze(rvEyeL(tt,2,:)), ...
                        'ZData',squeeze(rvEyeL(tt,3,:)));
                set(hRl,'XData',squeeze(rvEyeR(tt,1,:)), ...
                        'YData',squeeze(rvEyeR(tt,2,:)), ...
                        'ZData',squeeze(rvEyeR(tt,3,:)));
                count = 0;
                end

                if count > 5
                    set(hSl,'XData',NaN,'YData',NaN,'ZData',NaN);
                    set(hLl,'XData',NaN,'YData',NaN,'ZData',NaN);
                    set(hRl,'XData',NaN,'YData',NaN,'ZData',NaN);
                end
                 drawnow;
                 if recordMovie
                    writeVideo(v,getframe(gcf));
                 end
            end

            if recordMovie
             close(v);
            end
           
end

function rv = propagateAllStates(rv0,tau,mu)
         rv = NaN(numel(tau),6,size(rv0,1));
for ts=1:size(rv0,1)
     [~,rv(:,:,ts)] = pumpkyn.cr3bp.prop(tau,rv0(ts,:),mu);
end

end

function rv = propagateAllInitStates(rv0,tauVec,mu)
         rv = NaN(size(rv0));
for ts=1:size(rv0,1)
   [~,rvTmp] = pumpkyn.cr3bp.prop(tauVec,rv0(ts,:),mu);
    rv(ts,:) = rvTmp(end,:);
end

end