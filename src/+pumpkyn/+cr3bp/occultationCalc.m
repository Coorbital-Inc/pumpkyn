function [tauDur,pIdx,oIdx] = occultationCalc(tau,x,rP1,rP2,muStar)
%% Purpose:
%
%  This routine will determine the occultations of either primary body
%  in the CR3BP frame of reference. Be sure to include interpolated values for
%  the time and the states of the satellite as this routine detects rising
%  and falling edges.
%
%% Inputs:
%
%  tau                      [N x 1]                 Dimensionless Time
%
%  x                        [N x 6]                 Dimensionless Pos
%                                                   and Velocity of
%                                                   satellite in the CR3BP
%
%  rP1                      double                  Dimensionless Radius
%                                                   of Primary 1
%
%  rP2                      double                  Dimensionless Radius
%                                                   of Primary 2
%
%  muStar                   double                  Dimenionless Mass Ratio
%                                                   of both primaries
%                                                   muStar = m2/(m1+m2)
%
%% Outputs:
%
%  tauDur                   [M x 2]                 Intervals of
%                                                   dimensionless time
%                                                   where P(2|1) is blocking
%                                                   the LOS of the
%                                                   satellite to P(1|2)
%
%  pIdx                     [M x 1]                 Index of which
%                                                   primary body is
%                                                   blocking the satellite
%                                                   from seeing the other
%                                                   primary body.
%                                                   pIdx = 2, means P2 is
%                                                   blocking P1 from the
%                                                   SV.
%
%  oIdx                     [N x 1]                 Boolean Occultation
%                                                   index corresponding
%                                                   to an occultation
%                                                   occurance.
%                                                    true = occultation
%                                                   false = no occultation
%
%   
%
%% Revision History:
%  Darin C. Koblick                                         (c) 10-03-2025
%  Updated to include primary Idx                               06-10-2026
%  Copyright 2025 Coorbital, Inc.
%% --------------------- Begin Code Sequence ------------------------------
if nargin == 0
                        tau0 = 2*pi;
                          Np = 20;
                          pm = +1;
 [tau0, x0, mu, tStar, lStar] = pumpkyn.cr3bp.getTulip(tau0,Np,pm,1e-12);
                        %data = pumpkynPie.cr3bp.getCyclerData([],'3:4');
                        %  x0 = data(end).x;
                        %tau0 = data(end).tau;
                         %tau0 = 0.0482006363326616;
                         %  x0 = [ 0.996745053043675  0 0 0 1.15992877692432 0];
                     [tau,x] = pumpkyn.cr3bp.prop(tau0,x0,mu);
                         rP2 = 1737.1./lStar;
                         rP1 = 6378.0./lStar;
          [tauDur,pIdx,oIdx] = pumpkyn.cr3bp.occultationCalc(tau,x,rP1,rP2,mu);
                      occDur = diff(tauDur,1,2).*tStar./3600;  %Hrs
                      totOcc = sum(occDur);

               figure('color',[1 1 1]);
               plot3(x(:,1),x(:,2),x(:,3),'.k'); hold on;
               %plot3(x(oIdx_i,1),x(oIdx_i,2),x(oIdx_i,3),'.r');
               [xS,yS,zS] = sphere(30);
               surf(xS.*rP2 + 1 - mu, yS.*rP2, zS.*rP2,'faceColor','k');
               surf(xS.*rP1 - mu, yS.*rP1,zS.*rP1,'faceColor','b');
               %Interpolate the sv based on the tauDur:
                for td=1:size(tauDur,1)
                    x_o = interp1(tau,x,linspace(tauDur(td,1),tauDur(td,2),60)','spline');
                    plot3(x_o(:,1),x_o(:,2),x_o(:,3),'-r','linewidth',2);
                end
               axis equal;
               set(gca,'clipping','off');
                  return;
end

%% Compute Min Acceptable Arclength Based on Primary Radius:
    sMin = min([rP1,rP2]./4);

%% Interpolate SV based on arclength:
[sTot,s] = pumpkyn.util.arclength(x(:,1:3),2);

%% Form an interpolation vector w sufficient resolution to catch occultation:
  nArc  = max(2,ceil(sTot./sMin)+1);
    s_i = linspace(0,sTot,nArc)';
    x_i = interp1(s,x,s_i,'spline');
  tau_i = interp1(s,tau,s_i,'linear');

%% Vector Calculation:
    posP1 = [-muStar,0, 0];             % First Primary postion
    posP2 = [1-muStar,0,0];             % Second Primary position
    rS2P2 = posP2 - x_i(:,1:3);         % SV -> P2 Line of Sight
    rS2P1 = posP1 - x_i(:,1:3);         % SV -> P1 Line of Sight
 rS2P2Mag = pumpkyn.util.vmag(rS2P2,2);
 rS2P1Mag = pumpkyn.util.vmag(rS2P1,2);

%% Determine LOS obstruction:
  phi = pumpkyn.util.bsxAng(rS2P2,rS2P1,2);       % Angle between LOS vectors

%% Determine relative half-angles of P1 and P2:
phiP2 = asind(min(1,rP2./rS2P2Mag));  %Half-angle of P2 rel 2 SV
phiP1 = asind(min(1,rP1./rS2P1Mag));  %Half-angle of P1 rel 2 SV

%% Is the LOS occulted by either primary?
oIdx_i =  phi <= (phiP2 + phiP1);

%% Which primary is in front of the other?
                                pIdx_i = zeros(size(tau_i));
pIdx_i(oIdx_i & (rS2P2Mag < rS2P1Mag)) = 2;
pIdx_i(oIdx_i & (rS2P1Mag < rS2P2Mag)) = 1;

%% Determine arclength occultation intervals for each primary:
    dOcc = diff([false; pIdx_i == 2; false]);
idxStart = find(dOcc == 1);
 idxStop = find(dOcc == -1) - 1;
tauDurP2 = [tau_i(idxStart),tau_i(idxStop)]; %#ok<FNDSB>
    dOcc = diff([false; pIdx_i == 1; false]);
idxStart = find(dOcc == 1);
 idxStop = find(dOcc == -1) - 1;
tauDurP1 = [tau_i(idxStart),tau_i(idxStop)]; %#ok<FNDSB>
  tauDur = [tauDurP1; tauDurP2];

 %% Empty Case
if isempty(tauDur)
  tauDur = [0 0];
    pIdx = zeros(1,1);
    oIdx = false(size(tau));
    return;
end

 %% Sort based on midpoint time:
      tauMean = mean(tauDur,2);
[tauMean,idx] = sort(tauMean);
       tauDur = tauDur(idx,:);

%% Assign one blocker index per interval:
    pIdx = interp1(tau_i,pIdx_i,tauMean,'nearest');

%% Map occultation flag back to original tau samples:
    oIdx = logical(interp1(tau_i,double(oIdx_i),tau,'nearest'));
end
