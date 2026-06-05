function propIMSV2()

% Propagate the SV2 constellation (ELFOs) around the moon
% 1). convert from kepler to cr3bp
% 2). use cr3bp to propagate for a sidereal month

sma0  = [12002.3,12005.1,11988.1,12013.8,11997.3]';
inc0  = [  56.24,  56.14,  56.71,  55.84,  57.24]';
ecc0  = [  0.678,  0.706,  0.695, 0.683,   0.698]';
RAAN0 = [  40.10,  41.18, 161.41, 160.82,  279.93]'; 
argp0 = [  90.03,  89.64,  90.32,  90.69,   89.84]';
   M0 = [   0.62, 183.51,  60.37, 236.41,  114.83]';

%Convert from mean anomaly to true anomaly:
  nu0 = pumpkyn.util.mean2True(M0,ecc0);

%Constants:
muStar = 0.012150585609624;        %Mass ratio
     G = 6.67384e-20;              %Gravitational Constant
     M = 5.9736E24 + 7.35E22;      %Char mass (sum of primaries)
 lStar = 389703;
 tStar = 382981.289129055;

%Convert from keplerian to CR3BP:  
    oev0 = [sma0, ecc0, inc0.*pi/180, argp0.*pi/180, RAAN0.*pi/180, nu0.*pi/180];
      x0 = pumpkyn.cr3bp.fromOrb(zeros(size(oev0,1),1),oev0,M,muStar,tStar,lStar,2);
     tau = linspace(0,2*pi,3600)';
       x = NaN(numel(tau),6,size(x0,1));

%Propagate for a full month       
for ts=1:size(x0,1)     
    [~,x(:,:,ts)] = pumpkyn.cr3bp.prop(tau,x0(ts,:),muStar);
end

figure('color',[1 1 1]);     
plot3(squeeze(x(:,1,:)),squeeze(x(:,2,:)),squeeze(x(:,3,:)));
hold on;
set(gca,'ColorOrderIndex',1);
for ts=1:size(x,3)
    plot3(x(1,1,ts),x(1,2,ts),x(1,3,ts),'.','markersize',25);
end
axis equal;

end