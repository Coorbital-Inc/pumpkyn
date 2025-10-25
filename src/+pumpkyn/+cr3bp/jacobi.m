function [C,U] = jacobi(x0,muStar)
%% Purpose:
%
%  This routine will compute the jacobi constant for a given
%  periodic orbit or zero velocity surface.
%
%% Inputs:
%
%  x0                   [N x 6]             Dimensionless States (pos/vel)
%
%  muStar               double              Mass Ratio of Primaries
%                                           muStar = mu2/(mu1+mu2)
%
%% Outputs:
%
%  C                    [N x 1]             Jacobi Constant
%
%  U                    [N x 1]             Pseudo-Potential
%
%% Revision History:
%  Darin C. Koblick                                              09/15/2025
%  Copyright 2025 Coorbital, Inc.
%% --------------------------- Begin Code Sequence ------------------------
if nargin == 0
   %Test:
       x0 = [8.2339081983651485E-1	-1.9017764504099543E-28	9.8941366235910004E-4	...
            -2.3545391932685812E-15	1.2634272983881797E-1	2.2367029429442455E-16];
   muStar = 1.215058560962404E-2;
        C = pumpkyn.cr3bp.jacobi(x0,muStar);
        %3.1743435193301202E+0
        return;
end
 x = x0(:,1);
 y = x0(:,2);
 z = x0(:,3);
r1 = sqrt(       (x + muStar).^2 + y.^2 + z.^2 );
r2 = sqrt( (x + muStar - 1).^2 + y.^2 + z.^2 );
 U = 0.5.*(x.^2 + y.^2) + (1 - muStar)./r1 + muStar./r2;
V2 = sum(x0(:,4:6).^2,2);
 C = 2.*U - V2;
end