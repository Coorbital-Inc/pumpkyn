function gdopOpt()
% Purpose:
%
%  This is an example of how to optimize GDOP using nonlinear optimization
%  using multiple pattern searches with the Global Optimization Toolbox.
%
%% Revision History:
%  Darin C. Koblick                                             09-24-2025
%% ------------------------- Begin Code Sequence --------------------------

%% Constants:
useParallel = false;           %Serial avoids overhead for this fast objective
     tStar = 382981.289129055;
        Np = 7;                 %Number of petals
      dtau = 15*60/tStar;       %Time step size (nondimensional)
  minElAng = 0*pi/180;          %Min Elevation Angle (rad)
    llaObs = [-90*pi/180,0,0];    %Moon observer lat/lon/alt
      Nsat = 6;             %Total number of satellites
 
%% Set the GDOP objective:
     gdopLimit = 6;             %Max GDOP counted as successful coverage

%% Set the search budget:
       Nstart = 64;            %Number of random constellations to rank
      maxTime = 12;             %Total search time (sec), including ranking
 maxStartTime = 2;              %Maximum refinement time per start (sec)

%% Form the basis of a tulip-shaped orbit:
                    tau0 = (5/6)*2*pi;
                      pm = -1;                      
[tau,rTgt,~,mu,lStar] = pumpkyn.cr3bp.tulipConstellation(Np,tau0,1,pm,0);
                     rP2 = [1-mu,0,0];
                   rvObs = pumpkyn.cr3bp.fromLLA([],llaObs,mu,lStar,2,2);
                    rObs = rvObs(:,1:3);
                    
%Check that we have a good solution:
figure('color',[1 1 1]);
plot3(rTgt(:,1),rTgt(:,2),rTgt(:,3),'k'); hold on;
plot3(1-mu,0,0,'.k','markersize',15);
axis equal;
grid on;
                    
%% Cache the position spline once for all objective evaluations:
% Keep the true propagated period even when the sampling grid stops short.
           tau0 = tau(end);
          tau_i = (0:dtau:tau0)';
      tauSpline = unique([tau_i; tau0]);
           rTgt = interp1(tau,rTgt,tauSpline,'spline');
 positionSpline = spline(tauSpline,rTgt.');
            tau = tau_i;

%% Set the number of satellites and phases to optimize:
 fixFirstPhase = false;          %Fix satellite 1 at zero; optimize Nsat-1 phases
%Set fixFirstPhase=false to optimize all Nsat phases.
%For a fixed observer over one orbit, a common phase shift only shifts time.
%Discrete time sampling can still give slightly different GDOP percentages.
   fixedPhases = [];
   
if fixFirstPhase
   fixedPhases = 0;
end
          Nvar = Nsat-numel(fixedPhases);
            lb = zeros(Nvar,1);
            ub = ones(Nvar,1);
%Search phases in [0,1]; convert to orbit time for the objective:
     objective = @(x)computeGDOPObjective([fixedPhases,tau0.*x(:)'], ...
                         tau,tau0,positionSpline,rP2,rObs,minElAng,gdopLimit);
    p = gcp('nocreate');
   if isempty(p)
    useParallel = false;
   end

%% Rank random constellations:
  searchTimer = tic;
       xStart = rand(Nstart,Nvar);
       scores = Inf(Nstart,1);
for kk=1:Nstart
   scores(kk) = objective(xStart(kk,:));
   if toc(searchTimer) >= maxTime
       break;
   end
end
    [~,order] = sort(scores);
      bestObj = scores(order(1));
         xSol = xStart(order(1),:);

%% Refine the best starting points until the time budget is used:
         opts = optimoptions('patternsearch','Display','off', ...
                'UseParallel',useParallel,'MaxFunctionEvaluations',1200, ...
                'InitialMeshSize',0.15,'MeshTolerance',1e-4);
for kk=1:Nstart
    remaining = maxTime-toc(searchTimer);
    if remaining <= 0 || bestObj == 0
        break;
    end
 opts.MaxTime = min(maxStartTime,remaining);
     [x,fval] = patternsearch(objective,xStart(order(kk),:), ...
                             [],[],[],[],lb,ub,[],opts);
    if fval < bestObj
      bestObj = fval;
         xSol = x;
        fprintf(1,'Start %d: %0.2f%% coverage (%0.1f sec)\n', ...
                   kk,100.*(1-bestObj),toc(searchTimer));
    end
end

%Convert the best phases back to orbit time:
         xSol = tau0.*xSol;

%Include the fixed satellite in the final geometry and GDOP statistics:
pseudoAnomalyVec = sort([fixedPhases,xSol(:)']);
 [obj,rTgt,gdop] = computeGDOPObjective(pseudoAnomalyVec,tau,tau0,positionSpline,rP2,rObs,minElAng,gdopLimit);
 
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
ylim([0 gdopLimit]);
set(gca,'xTick',0:1:max(tau.*tStar./86400)+1);

%Output the GDOP Statistic:
fprintf(1,'%0.2f%% of time GDOP <= %g\n',100.*sum(gdop<=gdopLimit)./numel(gdop),gdopLimit)

end


function [obj,rTgt,gdop] = computeGDOPObjective(pseudoAnomalyVec,tau,tau0,positionSpline,rP2,rObs,minElAng,gdopLimit)
%Purpose:
%
% Minimize the fraction of time GDOP exceeds gdopLimit or is unavailable.

%Sort pseudoAnomalyVec and make sure it's unique:
pseudoAnomalyVec = unique(pseudoAnomalyVec);

rTgt = [];
gdop = NaN(numel(tau),1);

if numel(pseudoAnomalyVec) < 4
    obj = 1;
    return;
end

       tau_i = mod(tau + pseudoAnomalyVec(:)',tau0);
        rTgt = permute(reshape(ppval(positionSpline,tau_i(:).'), ...
                    3,numel(tau),numel(pseudoAnomalyVec)),[2 1 3]);
      visible = pumpkyn.cr3bp.elAng(rObs,rTgt,rP2,2,minElAng);
      maskIdx = permute(~visible,[1 3 2]);
          dop = pumpkyn.cr3bp.dop(rObs,rTgt,maskIdx,2);
         gdop = dop(:,1);
          obj = 1 - sum(gdop <= gdopLimit)./numel(gdop);
end
