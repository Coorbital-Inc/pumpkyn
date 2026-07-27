function [r,v] = moonPosVel(jd)
%% Purpose:
%
%  High precision Earth Moon position and velocity routine adapted 
%  from the Mooon algorithm outlined in vallado (alg. 31) pgs. 288
%  and also sumarized in "An Alternative Lunar Ephemeris Model 
%  for On-Board Flight Software Use" Review of Current Models. Pgs. 176-177
%  David G. Simpson, NASA Goddard Space Flight Center. 1999. RMS error in 
%  lunar position of 0.1 deg, with a maximum of 0.35 deg.
%
%  Output vectors are J2000 Mean Equator of Date (MOD) reference 
%  frame.  Velocity is computed from an analytical derivative of the 
%  position vector with respect to julian centuries, then converted to 
%  seconds before output.  This routine has been optimized for extreme 
%  computational speed to support mobile processing applications.
%
%  ECI Position vector output verified w Vallado and Astronomical Almanac
%  for April 28, 1994 00:00:00 UTC. Velocity vector output verified against 
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
%                                                   Vector of the moon 
%                                                   in km
%
%  v                            [N x 1]             MOD J2000 Velocity 
%                                                   Vector of the moon
%                                                   in km/s
%
%% Revision History:
%  Darin C. Koblick                                          (c) 07-27-26
%% --------------------------- Begin Code Sequence ------------------------
if nargin == 0
       jd = juliandate('04/28/1994 00:00:00');
       dt = 60*60;
        t = (0:dt:86400*27.5*5)';
    [r,v] = moonPosVel(jd+t./86400);
    figure('color',[1 1 1]);
    subplot(1,2,1);
    plot(t./86400,r,'-','linewidth',4);
    grid on;
    subplot(1,2,2);
    plot(t./86400,v,'-','linewidth',4); hold on;
    plot((t(1:end-1)+dt/2)./86400,diff(r,1)./dt,'-k','linewidth',2);  
    grid on;
    figure('color',[1 1 1]);
    plot3(r(:,1),r(:,2),r(:,3),'-k'); hold on;
    plot3(0,0,0,'.b','markersize',60);
    plot3(r(1,1),r(1,2),r(1,3),'.r');
    axis equal;
    grid off;
    axis off;
    
    %Determine the title angle:
    hMoon = cross(r,v,2);
    theta = bsxAng(hMoon,[0 0 1],2);
    
    figure('color',[1 1 1]);
    plot(t./86400,theta);
    grid on;
    ylabel('Incination [deg]');
    
    r = [];
    v = [];
    return;
end

%% Constant Parameters:
        ER2km = 6378.137;         %Equatorial Radius IAU 1976 Value: km/ER
         tc2s = 36525*24*3600;                                  %sec/C
        dg2rd = pi/180;                                         %rad/deg
      
%% Convert Julian Day to Julian Centuries       
         ttdb = (jd - 2451545.0)./36525.0;                      %C
       
%% Celestial Parameters:       
    eclplong = 218.32  + 481267.8813.*ttdb ...
                       + 6.29.*sind( 134.9 +477198.85.*ttdb) ...
                       - 1.27.*sind(  259.2 -413335.38.*ttdb) ...
                       + 0.66.*sind(  235.7 +890534.23.*ttdb) ...
                       + 0.21.*sind(  269.9 +954397.70.*ttdb) ...
                       - 0.19.*sind(  357.5 + 35999.05.*ttdb) ...
                       - 0.11.*sind(  186.6 +966404.05.*ttdb);   % deg
    cosEclLon = cosd(eclplong);    
    sinEclLon = sind(eclplong);
                
      eclplat =   5.13.*sind(93.3  +483202.03.*ttdb) ...
                + 0.28.*sind(228.2 +960400.87.*ttdb) ...
                - 0.28.*sind(318.3 +  6003.18.*ttdb) ...
                - 0.17.*sind(217.6 -407332.20 *ttdb);            % deg
            
    cosEclLat = cosd(eclplat);
    sinEclLat = sind(eclplat);
    
 %Horizontal Paralax:     
      hzparal =  0.9508  +  0.0518.*cosd(134.9 +477198.85.*ttdb) ...
                         +  0.0095.*cosd(259.2 -413335.38.*ttdb) ...
                         +  0.0078.*cosd(235.7 +890534.23.*ttdb) ...
                         +  0.0028.*cosd(269.9 +954397.70.*ttdb);% deg 
              
    obliquity = 23.439291  - 0.0130042.*ttdb;                    % deg               
       cosObl = cosd(obliquity);
       sinObl = sind(obliquity);
       
%% Determine the J2000 MOD Position Vector of the Moon (km):      
            l = cosEclLat.*cosd(eclplong );
            m = cosObl.*cosEclLat.*sinEclLon ...
                - sinObl.*sinEclLat;
            n = sinObl.*cosEclLat.*sinEclLon ...
                + cosObl.*sinEclLat;
         magr = 1./sind(hzparal);                                 %ER
            r = ER2km.*magr.*[l,m,n];                             %km                
if nargout == 2
%% Determine the J2000 MOD Velocity Vector of the Moon (km/s):            
   hzparalDot = -0.0518.*477198.85.*sind(134.9 +477198.85.*ttdb) + ...
                +0.0095.*413335.38.*sind(259.2 -413335.38.*ttdb) + ...
                -0.0078.*890534.23.*sind(235.7 +890534.23.*ttdb) + ...
                -0.0028.*954397.70.*sind(269.9 +954397.70.*ttdb); %deg/C
            
   hzparalDot = hzparalDot.*dg2rd.*dg2rd;                         %rad/C         
                
   eclplatDot = +5.13.*483202.03.*cosd(93.3  +483202.03.*ttdb) + ...
                +0.28.*960400.87.*cosd(228.2 +960400.87.*ttdb) + ...
                -0.28.*6003.18.*cosd(318.3 +  6003.18.*ttdb)   + ...
                +0.17.*407332.20.*cosd(217.6 -407332.20 *ttdb);   %deg/C
            
   eclplatDot = eclplatDot.*dg2rd.*dg2rd;                         %rad/C
                 
  eclplongDot = 481267.8813./dg2rd + ...         
              + 6.29.*477198.85.*cosd( 134.9 +477198.85.*ttdb) + ...
              + 1.27.*413335.38.*cosd(  259.2 -413335.38.*ttdb) + ...
              + 0.66.*890534.23.*cosd(  235.7 +890534.23.*ttdb) + ...
              + 0.21.*954397.70.*cosd(  269.9 +954397.70.*ttdb) + ...
              - 0.19.*35999.05.*cosd(   357.5 + 35999.05.*ttdb) + ...
              - 0.11.*966404.05.*cosd(  186.6 +966404.05.*ttdb);  %deg/C
          
  eclplongDot = eclplongDot.*dg2rd.*dg2rd;                        %rad/C
 obliquityDot = -0.0130042.*dg2rd;                                %rad/C        
          
      magrDot = -hzparalDot.*cosd(hzparal)./sind(hzparal).^2;     %ER/C
      
           vx = (-magr.*eclplatDot.*sinEclLat.*cosEclLon + ...
                 -magr.*cosEclLat.*eclplongDot.*sinEclLon + ...
                 +magrDot.*cosEclLat.*cosEclLon); 
                            
           vy = (-magr.*obliquityDot.*sinObl.*cosEclLat.*sinEclLon + ...
                 -magr.*cosObl.*eclplatDot.*sinEclLat.*sinEclLon + ...
                 +magr.*cosObl.*cosEclLat.*eclplongDot.*cosEclLon + ...
                 +magrDot.*cosObl.*cosEclLat.*sinEclLon - ...
                 +magr.*obliquityDot.*cosObl.*sinEclLat - ...
                 +magr.*sinObl.*eclplatDot.*cosEclLat - ...
                 +magrDot.*sinObl.*sinEclLat); 
                            
           vz = (+magrDot.*sinObl.*cosEclLat.*sinEclLon + ...
                 +magr.*obliquityDot.*cosObl.*cosEclLat.*sinEclLon + ...
                 -magr.*sinObl.*eclplatDot.*sinEclLat.*sinEclLon + ...
                 +magr.*sinObl.*cosEclLat.*eclplongDot.*cosEclLon + ...
                 +magrDot.*cosObl.*sinEclLat + ...
                 -magr.*obliquityDot.*sinObl.*sinEclLat + ...
                 +magr.*cosObl.*eclplatDot.*cosEclLat);                 
                            
           v = (ER2km./tc2s).*[vx,vy,vz];                          %km/s
end
end