function M = monodromy(x0,tau0,muStar)
%% Purpose:
%
%  This routine will compute the monodromy matrix assoicated with 
%  the inital dimensionless states, x0, and tau0 of a periodic orbit
%
%% Inputs:
%
%  tau0                 double              Dimensionless Period of Orbit
%
%  x0                   [1 x 6]             Dimensionless States (pos/vel)
%
%  muStar               double              Mass Ratio of Primaries
%                                           muStar = mu2/(mu1+mu2)
%
%% Ouputs:
%
%      M                [6 x 6]            Monodromy Matrix corresponding
%                                          to the STM @ t = tau0
%                                            
%% Revision History:
%  Darin C. Koblick                                              09/22/2025
%  Copyright 2025 Coorbital, Inc.
%% --------------------------- Begin Code Sequence ------------------------
if nargin == 0
   %two-dimensional lyaponuv orbit bout L2 
   tau0 = 3.4344774781239868;
     x0 = [1.1130439425606091E+0	
           -6.0030029039812299E-28	
           8.9831983586449705E-33	
           4.0854090008117549E-15	
           2.1082532553834954E-1	
           7.5485662781839097E-33]';
   muStar = 1.215058560962404E-2;
        M = pumpkyn.cr3bp.monodromy(x0,tau0,muStar);
        return;
end
       PHI0 = eye(6);
   [~,xPHI] = pumpkyn.cr3bp.prop(tau0,[x0(:); PHI0(:)],muStar);
          M = reshape(xPHI(end,7:end),[6 6]);
end