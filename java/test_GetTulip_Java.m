function test_GetTulip_Java()
%% Purpose:
%
%  Test wrapper for JAVA, getTulip, source code.
%

javaaddpath(which('GetTulip.jar'));

%Get Interpolated Orbit Data:
     Np = 7;
   tau0 = 2*pi;
     pm = -1;
   data = GetTulip.getTulip(tau0, Np, pm);

%% Extract Data and plot:
  tau0 = data(1);
    x0 = data(2:7);
muStar = data(8);
 tStar = data(9);
 lStar = data(10);
 [tau, x] = pumpkyn.cr3bp.prop(tau0,x0(:)',muStar);
 figure('Color',[1 1 1]);
 plot3(x(:,1).*lStar,x(:,2).*lStar,x(:,3).*lStar,'-','LineWidth',3);
 axis equal; set(gca,'clipping','off'); 

end