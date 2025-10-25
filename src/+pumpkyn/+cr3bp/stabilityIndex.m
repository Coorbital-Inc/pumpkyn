function k = stabilityIndex(x0,tau0,muStar)
%% Purpose:
%
%  This routine will compute the stability index of periodic orbit.
%
%% Inputs:
%
%  x0                   [1 x 6]             Dimensionless States (pos/vel)
%
%  tau0                 double              Dimensionless Period of Orbit
%
%  muStar               double              Mass Ratio of Primaries
%                                           muStar = mu2/(mu1+mu2)
%
%% Outputs:
%
%   k                   double              Stability Index
%
%% Revision History:
%  Darin C. Koblick                                              09/15/2025
%  Copyright 2025 Coorbital, Inc.
%% --------------------------- Begin Code Sequence ------------------------
if nargin == 0
       x0 = [8.2339081983651485E-1	-1.9017764504099543E-28	9.8941366235910004E-4	...
            -2.3545391932685812E-15	1.2634272983881797E-1	2.2367029429442455E-16];
     tau0 = 2.7430007981241529E+0;
   muStar = 1.215058560962404E-2;
        k = pumpkyn.cr3bp.stabilityIndex(x0,tau0,muStar);
   %1.1804065333857600E+3
   return;
end
     M = pumpkyn.cr3bp.monodromy(x0,tau0,muStar);
lambda = eig(M);
     k = 0.5.*(lambda + 1./lambda);
     k = max(abs(real(k)));
end