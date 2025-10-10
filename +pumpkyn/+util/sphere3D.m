function [x,y,z,xS,yS,zS] = sphere3D(N)
%% Purpose:
%
%  This routine will compute a sphere using parametric equations which
%  should space the points out equally over the spherical surface
%
%% Source:
%  How to generate equidistributed points on the surface of a sphere.
%  Markus Deserno. 2004.
%
%% Inputs:
%
%  N                        integer                      Number of Points
%
%
%% Outputs:
%
%  x                        [N x 1]                      x cartesian
%                                                        component
%
%  y                        [N x 1]                      y cartesian
%                                                        component
%
%  z                        [N x 1]                      z cartesian
%                                                        component
%
%  xS                       [M x P]                     x cartesian
%                                                       components needed
%                                                       for a surface plot
%
%  yS                       [M x P]                     y cartesian
%                                                       components needed
%                                                       for a surface plot
%
%  zS                       [M x P]                     z cartesian
%                                                       components needed
%                                                       for a surface plot
%
%% Revision History:
%  Darin C. Koblick                                     (c) 09/30/2025
%  Copyright 2025 Coorbital, Inc.
%% --------------------- Begin Code Sequence ------------------------------
if nargin == 0
   [x,y,z,xS,yS,zS] = pumpkyn.util.sphere3D(2000);
   figure('color',[1 1 1]);
   surf(xS,yS,zS,'faceAlpha',1,'edgeColor','none'); hold on;
   plot3(x(:),y(:),z(:),'.k');
   axis equal;
   return;
end

         r = 1;
         a = (4.*pi.*r.^2)./N;
      M_nu = round(pi./sqrt(a));
      d_nu = pi./M_nu;
     d_phi = a./d_nu;
         m = 0:M_nu-1;
        nu = pi.*(m+0.5)./M_nu;
     M_phi = round(2.*pi.*sin(nu)./d_phi);
         n = (0:max(M_phi))';
       phi = 2.*pi.*n./M_phi;
         x = sin(nu).*cos(phi);
         y = sin(nu).*sin(phi);
         z = repmat(cos(nu),[size(phi,1) 1]);
[xS,yS,zS] = deal(x,y,z);

     %Only preserve upto (0:M_phi(tm)-1)'
     for tm=1:numel(M_phi)
                   idx = ismember(1:size(x,1),1:M_phi(tm))';
            x(~idx,tm) = NaN;
            y(~idx,tm) = NaN;
            z(~idx,tm) = NaN;
     end

          x = x(:);
          y = y(:);
          z = z(:);
x(isnan(x)) = [];
y(isnan(y)) = [];
z(isnan(z)) = [];

end