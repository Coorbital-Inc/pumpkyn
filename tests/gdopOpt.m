function gdopOpt()
% Purpose:
%
%  This is an example of how to optimize GDOP using nonlinear optimization
%  using the genetic algorithm global optimization toolbox.
%
%% Revision History:
%  Darin C. Koblick                                             09-24-2025
%% ------------------------- Begin Code Sequence --------------------------

%% Constants:
useParallel = true;
     tStar = 382981.289129055;
        Np = 7;                  %Number of petals
      Ntau = Np*180;              %Number of steps per revolution
 
%% Form the basis of a tulip-shaped orbit:
                    tau0 = (5/6)*2*pi;
                      pm = -1;                      
[tau,rTgt,vTgt,mu,lStar] = pumpkyn.cr3bp.tulipConstellation(Np,tau0,1,pm,0);
                     rP2 = [1-mu,0,0];
                    rObs = rP2 + [0 0 -1738.1./lStar];
                    
%Check that we have a good solution:
figure('color',[1 1 1]);
plot3(rTgt(:,1),rTgt(:,2),rTgt(:,3),'k'); hold on;
plot3(1-mu,0,0,'.k','markersize',15);
axis equal;
grid on;
                    
%% Interpolate so we have evenly spaced timestpes:
 tau_i = linspace(0,tau0,Ntau)';
  rTgt = interp1(tau,rTgt,tau_i,'spline');
  vTgt = interp1(tau,vTgt,tau_i,'spline');
   tau = tau_i;
   
%% Given a pseudo-anomaly [0 to tau0] vector, optimize the trajectory:
 Nsat = 5;
   lb = zeros(Nsat,1);
   ub = repmat(tau0,[Nsat 1]);
    p = gcp('nocreate');
   if isempty(p)
    useParallel = false;
   end
%% USe Genetic Optimization   
  opts = optimoptions('ga','Display','iter','MaxTime',10*60,'UseParallel', ...
                      useParallel,'MaxStallGenerations',20,'PopulationSize',2^9,'CrossOverFraction',0.7);
  xSol = ga(@(x)computeGDOPObjective(x,tau,rTgt,vTgt,rP2,rObs),Nsat,[],[],[],[], ...
            lb,ub,[],opts);
       
%% Use PArticle Swarm:
%  opts = optimoptions('particleswarm','SwarmSize',200,'Display','iter','MaxTime',2*60,...
%         'useParallel',useParallel,'SelfAdjustmentWeight',1.2,'SocialAdjustmentWeight',1.2);
%  xSol = particleswarm(@(x)computeGDOPObjective(x,tau,rTgt,vTgt,rP2,rObs),Nsat, ...
%            lb,ub,opts);
       
      pseudoAnomalyVec = sort(xSol(:)');
 [obj,rTgt,vTgt,gdop] = computeGDOPObjective(pseudoAnomalyVec,tau,rTgt,vTgt,rP2,rObs);     
 
%Show Geometry: 
figure('color',[1 1 1]);
plot3(rTgt(:,1,1),rTgt(:,2,1),rTgt(:,3,1),'k'); hold on;
plot3(squeeze(rTgt(1,1,:)), ...
      squeeze(rTgt(1,2,:)), ...
      squeeze(rTgt(1,3,:)),'ok','markersize',8); 
axis equal;

%Show GDOP:
figure('color',[1 1 1]);
plot(tau.*tStar./86400,gdop,'k');
grid on;
xlabel('Time [Days]');
ylabel('GDOP');
ylim([0 6]);

set(gca,'xTick',0:1:max(tau.*tStar./86400)+1);
end


function [obj,rTgt,vTgt,gdop] = computeGDOPObjective(pseudoAnomalyVec,tau,rTgt,vTgt,rP2,rObs)
%Purpose:
%
% compute the percentage of time the GDOP is above 6 per revolution.
% Goal is to minimize this value as much as possible by re-arranging the
% pseudoAnomaly spacing of the satellites in the orbit.

%Sort pseudoAnomalyVec and make sure it's unique:
pseudoAnomalyVec = unique(pseudoAnomalyVec);

if numel(pseudoAnomalyVec) < 4
    obj = 1;
    return;
end

  [rTgt,vTgt] = genNewTrajAtPsuedoAnom(tau,rTgt,vTgt,pseudoAnomalyVec(:)');
          phi = pumpkyn.cr3bp.elAng(rObs,rTgt,rP2,2);
      maskIdx = permute(phi < 0*pi/180,[1 3 2]);  
          dop = pumpkyn.cr3bp.dop(rObs,rTgt,maskIdx,2);
         gdop = dop(:,1);
          obj = 1 - sum(gdop <= 6)./numel(gdop); %Percent GDOP above 6
end


function [r,v] = genNewTrajAtPsuedoAnom(tau,r,v,pseudoAnomaly)
          tau0 = tau(end);
         tau_i = tau + pseudoAnomaly;
             r = interp1(tau,r,mod(tau_i,tau0),'spline');
             v = interp1(tau,v,mod(tau_i,tau0),'spline');
             r = permute(r,[1 3 2]);
             v = permute(v,[1 3 2]);
end