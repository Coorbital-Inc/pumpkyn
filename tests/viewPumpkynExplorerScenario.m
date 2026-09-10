function result = viewPumpkynExplorerScenario(input,showPlots)
%% Purpose:
%
%  Import and evaluate a scenario exported by pumpkyn three-body explorer.
%  Use the pumpkyn CR3BP propagator, elevation-angle and DOP routines to
%  evaluate the constellation. Print summary statistics to the command
%  window and optionally display the constellation, GDOP and SV visibility.
%
%% Inputs:
%
%  input                char or string      Scenario JSON filename or its
%                                           text contents. Defaults to
%                                           'pumpkyn-scenario.json'.
%
%  showPlots            logical             Display figures (default true).
%                                           Set false for calculations only.
%
%% Outputs:
%
%  result               struct              Full-precision results:
%
%    scenario           struct              Imported parameters and states
%    tau                [Nt x 1]            Analysis times in TU
%    xObs               [1 x 6]             Stationary observer state, LU/TU
%    xSat               [Nt x 6 x Ns]       Phased spacecraft states, LU/TU
%    xOrbit             [Nt x 6 x N]        Unphased orbit states, LU/TU
%    dop                [Nt x 3]            GDOP, PDOP and TDOP
%    nLOS               [Nt x 1]            Number of spacecraft in LOS
%    stats              struct              gdopPercent, losPercent,
%                                           maxGapHours, meanLOSRangeKm
%    figures            graphics handles    Empty when showPlots is false
%
%  Uses the calculations in pumpkynExplorerScenario: periodic interpolation
%  of each saved orbit, sample-based coverage, and mean visible range. The
%  input analysis.compute_gdop flag does not skip this explicit evaluation.
%  The exported solver profile is retained in result.scenario; propagation
%  uses the existing settings in pumpkyn.cr3bp.prop.
%
%% Examples:
%
%  result = viewPumpkynExplorerScenario('pumpkyn-scenario.json');
%  result = viewPumpkynExplorerScenario(jsonText,false);
%
%% Revision History:
%  Copyright 2026 Coorbital, Inc.
%% -------------------------- Begin Code Sequence -------------------------
if nargin == 0
        viewPumpkynExplorerScenario('pumpkyn-scenario.json');
        result = [];
        return;
end

if nargin < 2
    showPlots = true;
end

%% Import scenario and constants:
     scenario = pumpkyn.util.importScenario(input);
       muStar = scenario.params.muStar;
        lStar = scenario.params.lStar;
        tStar = scenario.params.tStar;
           rM = scenario.params.rP2/lStar;
         xObs = scenario.xObs;
         tau0 = scenario.tau0;
       tauVec = (0:scenario.params.analysis.step:scenario.params.analysis.duration)';
            N = size(scenario.xSat,1);
gdopThreshold = scenario.params.analysis.gdop_threshold;
  requiredSVs = scenario.params.analysis.required_satellites;
     minElAng = deg2rad(scenario.params.observer.elevation_mask);

%% Propagate orbits and apply spacecraft offsets:
    xSat = cell(N,1);
    x3BO = NaN(numel(tauVec),6,N);
for k = 1:N
    [tauTmp,xTmp] = pumpkyn.cr3bp.prop(tau0(k),scenario.xSat(k,:),muStar);
          xSat{k} = interp1(tauTmp,xTmp, ...
                       mod(tauVec + scenario.dtau{k},tau0(k)),'spline');
      if size(xSat{k},3) > 1
           xSat{k} = permute(xSat{k},[1 3 2]);
      end
      x3BO(:,:,k) = interp1(tauTmp,xTmp,mod(tauVec,tau0(k)),'spline');
end
    xSat = cat(3,xSat{:});

%% Compute visibility and GDOP:
          phi = pumpkyn.cr3bp.elAng(xObs(:,1:3),xSat(:,1:3,:),[1-muStar,0,0],2);
          dop = pumpkyn.cr3bp.dop(xObs(:,1:3),xSat(:,1:3,:),permute(phi < minElAng,[1 3 2]),2);
         gdop = dop(:,1);
      gdopMed = median(gdop(isnumeric(gdop)));
 gdopFraction = mean(gdop < gdopThreshold);
  NsatsInView = sum(phi >= minElAng,3);
     LOSAvail = sum(NsatsInView >= requiredSVs,1)./numel(tauVec);
       maxGap = max(diff(unique([0; tauVec(gdop < gdopThreshold); tauVec(end)]),1,1))  ...
                - scenario.params.analysis.step;
       rng2O = pumpkyn.util.vmag(xSat(:,1:3,:) - xObs(:,1:3),2);
      rngLOS = rng2O(phi >= minElAng);
  meanLOSrng = mean(rngLOS);

%% Return full-precision results and print the summary:
            result.scenario = scenario;
                 result.tau = tauVec;
                result.xObs = xObs;
                result.xSat = xSat;
              result.xOrbit = x3BO;
                 result.dop = dop;
                result.nLOS = NsatsInView;
   result.stats.gdopPercent = gdopFraction*100;
    result.stats.losPercent = LOSAvail*100;
   result.stats.maxGapHours = maxGap*tStar/3600;
result.stats.meanLOSRangeKm = meanLOSrng*lStar;
             result.figures = gobjects(0);
fprintf('\nScenario: %s - %s | %d orbits | %d SVs\n', ...
    scenario.params.P1,scenario.params.P2,N,size(xSat,3));
fprintf('  %-18s: %s %%\n',['GDOP < ',num2str(gdopThreshold)], ...
    statText(result.stats.gdopPercent));
fprintf('  %-18s: %s %%\n',sprintf('LOS (%d+ SVs)',requiredSVs),statText(result.stats.losPercent));
fprintf('  %-18s: %s hrs\n','Max gap',statText(result.stats.maxGapHours));
fprintf('  %-18s: %s km\n\n','Mean LOS range',statText(result.stats.meanLOSRangeKm));
if ~showPlots
    return;
end

%% View the constellation:
[xM,yM,zM] = sphere(15);
result.figures(1) = figure('color',[1 1 1]);
surf(xM*rM + (1-muStar),yM*rM,zM*rM, ...
    'FaceColor','none','edgeColor','k','HandleVisibility','off'); hold on;
plot3(squeeze(x3BO(:,1,:)),squeeze(x3BO(:,2,:)),squeeze(x3BO(:,3,:)));
plot3(xObs(:,1),xObs(:,2),xObs(:,3),'.r','markersize',14);
plot3(squeeze(xSat(1,1,:)),squeeze(xSat(1,2,:)),squeeze(xSat(1,3,:)), ...
    '.k','markersize',15);
grid on; axis equal;
sNames = cellfun(@(o) o.name,scenario.orbits(:),'UniformOutput',false);
legend([sNames; 'Observer'; sprintf('SV 1-%d',size(xSat,3))],'Interpreter','none');
set(gca,'clipping','off');

%% View GDOP and satellites in LOS over time:
result.figures(2) = figure('color',[1 1 1]);
tiledlayout(2,1,'TileSpacing','compact','Padding','compact');
timeDays = tauVec*tStar/86400;
  axGDOP = nexttile;
plot(timeDays,gdop,'k','linewidth',2);
grid on; ylim([max(gdopMed-8,0),gdopMed+8]);
ylabel('GDOP');
   axLOS = nexttile;
stairs(timeDays,NsatsInView,'k','linewidth',1.5);
grid on; ylabel('SVs in LOS'); xlabel('Time [Days]');
yticks(min(NsatsInView):max(NsatsInView));
linkaxes([axGDOP,axLOS],'x');
end

function text = statText(value)
%Print three significant figures without scientific notation:
if isfinite(value) && value ~= 0
    text = sprintf('%.*f',max(0,2-floor(log10(abs(value)))),round(value,3,'significant'));
else
    text = num2str(value);
end
end
