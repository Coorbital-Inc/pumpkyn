
# Low Thrust Transfer from L1 Halo to a Pumpkin Orbit

This tutorial demonstrates how to use tfMin, a minimum time low\-thrust solver to transfer from an L1 halo orbit to a pumpkin orbit

# Characteristic mass (sum of primaries)
```matlab
                          M = 5.9736E24 + 7.35E22;  
```
# Initial Halo Orbit in L1:
```matlab
                         x0 = [8.2338739411138306E-1    0 6.9135470535882995E-3    0    1.2711946779802655E-1    0];
                       tau0 = 2.7433219717780992E+0;
                     muStar = 0.012150585609624;
                      lStar = 389703.264829278;
                      tStar = 382981.289129055;
             [tauInt,rvInt] = pumpkyn.cr3bp.prop(tau0,x0,muStar);
```
# Desired orbit to transfer to from interpolating seeds:
```matlab
                      Np = 7;
                      pm = +1;     
                    tau0 = 2*pi;
                  [~,x0] = pumpkyn.cr3bp.getTulip(tau0,Np,pm,1e-8);
          [tauTgt,rvTgt] = pumpkyn.cr3bp.prop(tau0,x0,muStar);
```
# Compute Lagrange points for reference:
```matlab
         Lpts = pumpkyn.cr3bp.lagrangePts(muStar);          
```
# Configure the Spacecraft Mass and Propulsion Parameters:
```matlab
       m0 = 1500;
       g0 = 9.80665 * tStar^2/(1000*lStar);        %m/s^2 -> ND
     Tmax = 4.5; %4.25;                                   %Newtons
```

```matlab
     Tmax = (Tmax/m0)*tStar^2/(lStar*1000);        %Max Thrust ND 
      Isp = 3000/tStar;    %Specific Impulse s -> ND
        c = Isp*g0;         %Exhuast Velocity (ND)
```
# Choose Boundary values:
```matlab
       rv0 = rvInt(52,:);
       rvf = rvTgt(160,:);  %160, %165
lambda0_tf = [  -2.43889808721205
          3.07944940888076
         -1.36965975864007
        -0.450614552645554
       -0.0443405922325243
         -0.38634620190255
1.24538681592255
2. 44736431592149];
```
# Solve the Problem w/ Minimum Transfer Time:
```matlab
sol_lambda0_tf = pumpkyn.cr3bp.tfMin(rv0,rvf,lambda0_tf,Tmax,c,muStar);
```
# Propagate Trajectory:
```matlab
 [tau,rv] = pumpkyn.cr3bp.tfMinProp(sol_lambda0_tf(8), ...
                  [rv0,1,sol_lambda0_tf(1:7)'],Tmax,c,muStar);
```
# Compute the Propellant Mass through the transfer:
```matlab
figure('color',[1 1 1]);
 mTot = rv(:,7).*m0;
mProp = m0 - mTot;
dVtot = c.*log(m0./mTot).*lStar./tStar;
subplot(1,2,1);
plot(tau.*tStar./86400,mProp);
xlabel('Transfer Time [Days]');
ylabel('Req Prop Mass [kg]');
grid on;
subplot(1,2,2);
plot(tau.*tStar./86400,dVtot);
xlabel('Transfer Time [Days]');
ylabel('\DeltaV_{tot} [km/s]');
grid on;
```
# Show the Trajectories
```matlab
     hIn = figure('color',[0 0 0]);
     pumpkyn.cr3bp.showMoon(lStar,muStar,hIn);
     plot3(Lpts(1,1),Lpts(1,2),Lpts(1,3),'.w','markersize',15); hold on;
     plot3(rvInt(:,1),rvInt(:,2),rvInt(:,3),'g');
     plot3(rvTgt(:,1),rvTgt(:,2),rvTgt(:,3),'r');
     plot3(rv(:,1),rv(:,2),rv(:,3),'w');
     axis equal;
     grid on; set(gca,'Clipping','off');
```
