function [r,v,a] = sunPosVel(jd)
%% Purpose:
%
%  Medium precision sun position and velocity routine adapted 
%  from the sun algorithm outlined in vallado (alg. 29) pgs. 267-280
%  Valid between 1950 - 2050.  Accuracy of apparent  coordinates is ~0.01 
%  degrees. Output vectors are J2000 Mean Equator of Date (MOD) reference 
%  frame.  Velocity is computed from an analytical derivative of the 
%  position vector with respect to julian centuries, then converted to 
%  seconds before output.  This routine has been optimized for extreme 
%  computational speed to support mobile processing applications.
%
%  ECI Position vector output verified w Vallado and Astronomical Almanac
%  for April 2, 2006 00:00:00 UTC. Velocity vector output verified against 
%  numerical derivative of position vector.
%
%
%% Inputs:
%
%  jd                           [N x 1]             Julian Day Vector
%
%
%% Outputs:
%
%  r                            [N x 1]             MOD J2000 Position
%                                                   Vector of the sun in km
%
%  v                            [N x 1]             MOD J2000 Velocity 
%                                                   Vector of the sun
%                                                   in km/s
%
%  a                            [N x 1]             MOD J2000 Acceleration
%                                                   Vector of the sun in
%                                                   km/s^2
%
%% Revision History:
%  Darin C. Koblick                                          (c) 06-26-2020
%  Darin C. Koblick    Added Acceleration Vector Output          01-22-2021
%% --------------------------- Begin Code Sequence ------------------------
if nargin == 0
       jd = pumpkyn.util.juliandate('04/02/1994 00:00:00');
        t = (0:60*60:86400*365*5)';
    tic;
    [r,v,a] = pumpkyn.util.sunPosVel(jd+t./86400);
    disp(toc);
    figure('color',[1 1 1]);
     subplot(1,3,1);
     plot(t./(86400*365.25),r,'-','linewidth',3); hold on;
     grid on;
     xlabel('Time [YRs]');
     ylabel('ECI J2000 Sun Position [km]');
     subplot(1,3,2);
     plot(t./(86400*365.25),v,'-','linewidth',3); hold on;
     plot(t(1:end-1)./(86400*365.25),diff(r,1,1)./diff(t,1),'-k');
     grid on;
     xlabel('Time [YRs]');
     ylabel('ECI J2000 Sun Velocity [km/s]');
     subplot(1,3,3);
     plot(t./(86400*365.25),a,'-','linewidth',3); hold on;
     plot(t(1:end-1)./(86400*365.25),diff(v,1,1)./diff(t,1),'-k');
     grid on;
     xlabel('Time [YRs]');
     ylabel('ECI J2000 Sun Acceleration [km/s^2]');
     %Show Earth Orbiting the Sun:
     figure('color',[1 1 1]);
     plot3(0,0,0,'.y','markersize',60); hold on;
     plot3(-r(:,1),-r(:,2),-r(:,3),'-k','linewidth',2);
     plot3(-r(1,1),-r(1,2),-r(1,3),'.b','markersize',15);
     axis equal;
     axis off;
     r = [];
     v = [];
     %Note: Astronomical Almanac lists the sun ICRS position vector for
     % April 2, 2006 00:00:00 UTC as [146,259,933 28595,947 12,397,430] km
    return;
end

%% Constant Parameters:
        AU2km = 149597870.0;                                       %AU/km
         tc2s = 36525*24*3600;                                     %C/sec
        dg2rd = pi/180;                                            %rad/deg
      
%% Convert Julian Day to Julian Centuries       
         tut1 = (jd - 2451545.0)./36525.0;                         %CC
       
%% Celestial Parameters:       
     meanlong = 280.460  + 36000.77.*tut1;                         %deg
         ttdb = tut1;
  meananomaly = 357.5277233  + 35999.05034.*ttdb;                  %deg
        cosMA = cosd(meananomaly);
        sinMA = sind(meananomaly);
     eclplong = meanlong + 1.914666471.*sinMA ...
              + 0.019994643.*sind(2.0.*meananomaly);               %deg
    cosEclLon = cosd(eclplong); 
    sinEclLon = sind(eclplong); 
    obliquity = 23.439291-0.0130042.*ttdb;                         %deg
       cosObl = cosd(obliquity);
       sinObl = sind(obliquity);
         magr = 1.000140612  - 0.016708617.*cosMA ...
              - 0.000139589.*cosd(2.0.*meananomaly);               %AU

%% Determine the J2000 MOD Position Vector of the Sun (km):                             
            r = AU2km.*magr.*[cosEclLon, ...
                              cosObl.*sinEclLon, ...
                              sinObl.*sinEclLon];                  %km
if nargout > 1                                           
%% Determine the J2000 MOD Velocity Vector of the Sun (km/s)
%Note: verified analytical calculations of derivative w/ numerical
%differentiation.
meananomalyDot =  35999.05034*dg2rd;                               %rad/C
   meanlongDot =  36000.77*dg2rd;                                  %rad/C
  obliquityDot = -0.0130042.*dg2rd;                                %rad/C
       magrDot = 0.016708617.*meananomalyDot.*sinMA + ...          %AU/C
                 2.*meananomalyDot.*0.000139589.*sind(2.0.*meananomaly);
   eclplongDot = meanlongDot + ...                                 %rad/C
                 1.914666471.*dg2rd.*meananomalyDot.*cosMA + ...
                 0.019994643.*dg2rd.*2.0.*meananomalyDot.*cosd(2.0.*meananomaly);        
            vx = -magr.*eclplongDot.*sinEclLon + ...
                  magrDot.*cosEclLon;                              %AU/C
            vy = -magr.*obliquityDot.*sinObl.*sinEclLon +  ...
                  magr.*cosObl.*eclplongDot.*cosEclLon + ...
                  magrDot.*cosObl.*sinEclLon;                      %AU/C
            vz = magr.*obliquityDot.*cosObl.*sinEclLon + ...
                 magr.*sinObl.*eclplongDot.*cosEclLon + ...
                 magrDot.*sinObl.*sinEclLon;                       %AU/C    
             v = (AU2km./tc2s).*[vx,vy,vz];                        %km/s
end

if nargout > 2
%% Determine the J2000 MOD Acceleration Vector of the Sun (km/s^2)
%Note: verified analytical calculations of second derivative w/ numerical
%differentiation.  
 eclplongDotDot =-1.914666471.*dg2rd.*meananomalyDot.*sinMA.*meananomalyDot - ...
                  0.019994643.*dg2rd.*2.0.*meananomalyDot.*sind(2.*meananomaly).*(2.*meananomalyDot);
     magrDotDot = 0.016708617.*meananomalyDot.*cosMA.*meananomalyDot + ...          
                  2.*meananomalyDot.*0.000139589.*cosd(2.*meananomaly).*(2.*meananomalyDot);
             ax = -magr.*eclplongDot.*cosEclLon.*eclplongDot + ...
                  -magr.*eclplongDotDot.*sinEclLon + ...
                  -magrDot.*eclplongDot.*sinEclLon + ...
                  -magrDot.*sinEclLon.*eclplongDot + ...
                   magrDotDot.*cosEclLon;
             ay = -magr.*obliquityDot.*sinObl.*cosEclLon.*eclplongDot +  ...
                  -magr.*obliquityDot.*cosObl.*obliquityDot.*sinEclLon +  ...
                  -magrDot.*obliquityDot.*sinObl.*sinEclLon + ...
                  -magr.*cosObl.*eclplongDot.*sinEclLon.*eclplongDot + ...
                   magr.*cosObl.*eclplongDotDot.*cosEclLon + ...
                  -magr.*sinObl.*obliquityDot.*eclplongDot.*cosEclLon + ...
                   magrDot.*cosObl.*eclplongDot.*cosEclLon + ...
                   magrDot.*cosObl.*cosEclLon.*eclplongDot + ...
                  -magrDot.*sinObl.*obliquityDot.*sinEclLon + ...
                   magrDotDot.*cosObl.*sinEclLon;
             az =  magr.*obliquityDot.*cosObl.*cosEclLon.*eclplongDot + ...   
                  -magr.*obliquityDot.*sinObl.*obliquityDot.*sinEclLon + ...
                   magrDot.*obliquityDot.*cosObl.*sinEclLon + ...
                  -magr.*sinObl.*eclplongDot.*sinEclLon.*eclplongDot + ...
                   magr.*sinObl.*eclplongDotDot.*cosEclLon + ...
                   magr.*cosObl.*obliquityDot.*eclplongDot.*cosEclLon + ...
                   magrDot.*sinObl.*eclplongDot.*cosEclLon + ...
                   magrDot.*sinObl.*cosEclLon.*eclplongDot + ...
                   magrDot.*cosObl.*obliquityDot.*sinEclLon + ...
                   magrDotDot.*sinObl.*sinEclLon;        
              a = (AU2km./tc2s.^2).*[ax,ay,az];                    %km/s^2
end

end