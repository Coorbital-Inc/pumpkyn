function gdop()

% This example demonstrates how we can compute the GDOP between a target
% and an observer in the cr3bp
                   tStar = 382981.289129055;
                      dt = 30*60./tStar;  %timestep
                    Npts = 500;
                    tau0 = (5/6)*2*pi;
                      pm = -1;
                      Np = 7;
                      
                    %6 Satellite Const   (0.0)
                     dtau =  [0.171755901906622
                             0.538936510214872
                             1.91530114835572
                             2.28537063968867
                             3.66203724308906
                             4.02638956895957]';
                        
%                   %5 Satellite Const (0.2794)    
%                   dtau = [1.49348712139373
%                           3.2267168526248
%                           3.57416662739935
%                           4.62537897604226
%                           4.98144670896738]';    
                   
%Propagate a set of tulip-shaped orbits:         
[tau,rTgt,vTgt,mu,lStar] = pumpkyn.cr3bp.tulipConstellation(Np,tau0,1,pm,dtau);
                      rP = 1738.1./lStar;   %Radius of P2
                      rE = 6378.137./lStar; %Radius of P1
                     rP2 = [1-mu,0,0];    %Location of P2 rel to barycenter
                   tau_i = linspace(0,tau(end),round(tau(end)/dt))';
                    rTgt = interp1(tau,rTgt,tau_i,'spline');
                    vTgt = interp1(tau,vTgt,tau_i,'spline');
                     tau = tau_i;
                     
%Compute the observers in the cr3bp:
     [rObs,lla,xyz] = pumpkyn.cr3bp.pointSphere(Npts,rP2,rP);
                lat = lla(:,1);
               rObs = permute(rObs,[4 2 3 1]); %[1 x 3 x 1 x Npts]
                phi = permute(pumpkyn.cr3bp.elAng(rObs,rTgt,rP2,2),[1 3 4 2]);   %[tau x Ntgts x Npts]
                idx = phi > 0*pi/180;  
%For each time step, determine the number of targets in view for each pt:
         NumVisTgts = permute(sum(idx,2),[1 3 2]); %[tau x Npts]
%For each time step, determine the DOP:         
                dop = pumpkyn.cr3bp.dop(rObs,rTgt,~idx,2);
                   
%Analyze lunar occultation:           
       rP1 = [-mu,0,0];
rTgt2Earth = rP1-rTgt;
 rTgt2Moon = rP2-rTgt;
     theta = squeeze(pumpkyn.util.bsxAng(rTgt2Earth,rTgt2Moon,2));   
   thetaP2 = atand(rP./pumpkyn.util.vmag(rTgt2Moon,2));
   thetaP1 = atand(rE./pumpkyn.util.vmag(rTgt2Earth,2));
    minAng = squeeze(thetaP2 + thetaP1);
    
    
%Look at south pole gdop region:
[val,idx] = min(lat);
gdop = squeeze(dop(:,1,idx(:)));
gdop = gdop(:,1);


%Determine the gap times between periods when GDOP > 6
idxGoodGDOP = gdop <= 6;
%   gapTimes = tau(idxGoodGDOP);
% GDOPOutage = diff(gapTimes);
%  idxOutage = GDOPOutage > dt + 1/tStar;

%Determine each outage interval:
            d = diff([NaN; idxGoodGDOP; NaN]);
    changeIdx = find(d ~= 0);
       starts = changeIdx(1:end-1);
         ends = changeIdx(2:end)-1;
       values = idxGoodGDOP(starts);
   outage_idx = [starts(values==0),ends(values==0)];
     gdop_idx = [starts(values==1),ends(values==1)];
outage_int_hr = (dt + tau(outage_idx(:,2)) - tau(outage_idx(:,1))).*tStar./3600;
  gdop_int_hr = (dt + tau(gdop_idx(:,2)) - tau(gdop_idx(:,1))).*tStar./3600;

  
%figure('color',[1 1 1]);


figure('color',[0 0 0]);
subplot(4,1,1);
plot(tau.*tStar./86400,gdop,'w'); hold on;
plot([tau(1),tau(end)].*tStar./86400,[6 6],'-','linewidth',2,'color',[1 0 0 0.5]);
ylabel('GDOP');
grid on;
ylim([0 7]);
title('\color{white}GDOP At South Pole');
set(gca,'color','k','xcolor','w','ycolor','w');

subplot(4,1,2); hold on;
for tt=1:size(gdop_idx,1)
    width = (tau(gdop_idx(tt,2)) - tau(gdop_idx(tt,1))).*tStar./86400;
    rectangle('Position', [tau(gdop_idx(tt,1)).*tStar./86400, 0, width, gdop_int_hr(tt)], ...
              'FaceColor', [0 1 0], 'EdgeColor', 'none');
end

for tt=1:size(outage_idx,1)
        width = (tau(outage_idx(tt,2)) - tau(outage_idx(tt,1))).*tStar./86400;
       rectangle('Position', [tau(outage_idx(tt,1)).*tStar./86400, 0, width, outage_int_hr(tt)], ...
                 'FaceColor', [1 0 0], 'EdgeColor', 'none'); 
end
grid on;
ylabel('Duration [Hrs]');
title(['\color{white}Max Outage = ',sprintf('%2.1f',max(outage_int_hr)),' Hrs']);
set(gca,'color','k','xcolor','w','ycolor','w');

%Look at number of visible satellites:
subplot(4,1,3);
Nvis = NumVisTgts(:,idx(:));
plot(tau.*tStar./86400,Nvis(:,1),'w');
ylabel('N_{sat}');
grid on;
ylim([numel(dtau)-2 numel(dtau)]);
title('\color{white}Number of Visible Satellites');
set(gca,'color','k','xcolor','w','ycolor','w','yTick',(-2:0)+numel(dtau));


%Look at lunar occlusions:
subplot(4,1,4);
plot(tau.*tStar./86400,min(theta,[],2),'w'); hold on;
plot(tau.*tStar./86400,max(minAng,[],2),'r');
legend('\color{white}Min Earth-Moon Angle','\color{red}Occultation Angle','color','none');
xlabel('Time [Days]');
ylabel('Earth-Moon Angle [deg]');
grid on;
title('\color{white}Lunar Occultation');
set(gca,'color','k','xcolor','w','ycolor','w');



                
%Animate the sequence w/ GDOP Overlay
            h = pumpkyn.cr3bp.showEarth(lStar,mu);
            pumpkyn.cr3bp.showMoon(lStar,mu,h);
            %hM = moon3D(rP2,true,1./lStar);
%thisNumVisTgts = reshape(NumVisTgts(1,:),size(x));
      tickVals = [0,4,6,8,20];
          gdop = reshape(dop(1,1,:),size(xyz(:,:,1)));
         % Use the Moon's axes explicitly
         ax = gca;
         hObs = surf(ax,(xyz(:,:,1) - rP2(1)).*1.2 + rP2(1), ...
                         xyz(:,:,2).*1.2, ...
                         xyz(:,:,3).*1.2, ...
                         gdop,'faceAlpha',1,'edgeColor','none'); 
         colormap(ax,jet(numel(tickVals)));
         %shading(ax,'interp');
          hColor = colorbar(ax); 
          caxis(ax, [tickVals(1), tickVals(end)]);
          hColor.Limits = [min(tickVals) max(tickVals)];
          hColor.Ticks = tickVals;
          hColor.Color = 'w';
          hColor.Label.String = 'gdop';
          for tt=1:size(rTgt,3)
            plot3(rTgt(:,1,tt),rTgt(:,2,tt),rTgt(:,3,tt),'w');
          end
          hTgts = plot3(squeeze(rTgt(1,1,:)),squeeze(rTgt(1,2,:)),squeeze(rTgt(1,3,:)),'.w','markersize',20);
          axis equal off;
          set(gca,'color','k','clipping','off');
          view(90,0);
          xlim([-1 +1].*max(max(abs(squeeze(rTgt(:,1,:))))));
          ylim([-1 +1].*max(max(abs(squeeze(rTgt(:,2,:))))));
          zlim([-1 +1].*max(max(abs(squeeze(rTgt(:,3,:))))));
          
for tt=1:numel(tau)      
     %thisNumVisTgts = reshape(NumVisTgts(tt,:),size(x));
     gdop = reshape(dop(tt,1,:),size(xyz(:,:,1)));
     set(hObs,'cData',gdop);
     set(hTgts,'xData',squeeze(rTgt(tt,1,:)), ...
               'yData',squeeze(rTgt(tt,2,:)), ...
               'zData',squeeze(rTgt(tt,3,:)));
     drawnow;
end


end