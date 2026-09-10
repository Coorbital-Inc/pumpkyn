function prop3BodyTraj()
%% Purpose:
%
%  This routine will take a trajectory from another application
%  and propagate it using the pumpkyn toolkit.
%

%Constants:
    muStar = 0.012150585609624;
     lStar = 389703.264829278;
     tStar = 382981.289129055;
        rM = 1737.1;
      tau0 = 24.735167667217148.*86400;
        x0 = [397910.7958059001	0	-22205.85237848581	0	0.16543480811895828	0];
      [tau0,x0] = pumpkyn.cr3bp.nondim(tau0,x0,tStar,lStar,2);

%Compute the statistics:
data = pumpkyn.cr3bp.orbitProperties(x0,tau0,muStar,lStar);

%Print out values:
fprintf(1,'Jacobi = %f [km^2/s^2]\n',data.Jacobi.*lStar.^2./tStar.^2);
fprintf(1,'Stability Index = %f\n',data.StabilityIndex);
fprintf(1,'Perilune = %f [km]\n',  data.Perilune.*lStar - rM);
fprintf(1,'Apolune = %f [km]\n',   data.Apolune.*lStar - rM);
fprintf(1,'Max Lunar Occultation = %f [days]\n',data.MaxLunarOcc.*tStar/86400);
fprintf(1,'Total Lunar Occultation = %f [days]\n',data.TotLunarOcc.*tStar/86400);

% figure('color',[1 1 1]);
% plot3(data.x(:,1),data.x(:,2),data.x(:,3));
% grid on;
% axis equal;

end