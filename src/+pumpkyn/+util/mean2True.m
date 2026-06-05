function TA = mean2True(MA,e)
%% Purpose:
%  Convert mean anomaly to true anomaly
%
%% Inputs:
%  MA                   [N x 1]                          Mean Anomaly (deg)
%
%  e                    [N x 1]                          Eccentricity of
%                                                        Orbit
%
%% Ouputs:
%  TA                   [N x 1]                          True Anomaly (deg)
%
%% Revision History:
%  Darin C. Koblick                                          (c) 08/14/2014
%% ----------------- Begin Code Sequence ----------------------------------
if nargin == 0
   MA = (0:360)';
    e = ones(size(MA)).*0.1;
    tic;
    TA = pumpkyn.util.mean2True(MA,e);
    toc
    return;
end
tol = 1e-13;
%% Step One: Find the Eccentric Anomaly (E) from the Mean Anomaly:
     MA = MA.*(pi/180);

%Vectorized newton solver approach:
     computeSol = @(E)E-e.*sin(E)-MA;
computeSolDeriv = @(E)1-e.*cos(E);
  computeUpdate = @(E)E-computeSol(E)./computeSolDeriv(E);
              E = MA;               %Initial Guess
         E_next = computeUpdate(E); %Initial Update to Guess
while (any(abs(E_next-E) > tol)) 
        E = E_next;
   E_next = computeUpdate(E);        
end
E = E_next;

% %Fzero newton solver approach:     
% options = optimset(@fzero);
%  findEA = @(e,MA)fzero(@(x)x-e.*sin(x)-MA,MA,options);
%       E = arrayfun(findEA,e,MA);

%% Step Two: Find the True Anomaly: 
     TA = Inf(size(MA));
    idx = e <= 1;
if any(idx)
    TA(idx) = 2.*atan2(sqrt(1+e(idx)).*sin(E(idx)./2),sqrt(1-e(idx)).*cos(E(idx)./2)).*(180/pi);
end

end