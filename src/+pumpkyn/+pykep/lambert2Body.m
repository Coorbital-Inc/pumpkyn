function [v0,vf] = lambert2Body(r0,rf,dt,mu,dir,nmax,dim3)
%% Purpose:
%
%  This is a wrapper interface for the pykep implementation of the
%  multi-revolution lambert solver. Note that this is applicable to
%  two-body equations of motion (keplerian).
%
%% Inputs:
%
%  r0                       [1 x 3]             Initial Position Vector
%                                               at t = 0
%
%  rf                       [1 x 3]             Final Position Vector
%                                               at t = dt
%
%  dt                       double              Time of flight (s)
%
%  mu                       double              Standard Gravitational
%                                               Parameter
%
%  dir                      integer             +1 = prograde transfer
%                                               -1 = retrograde transfer
%
%  nmax                     integer             Maximum number of
%                                               revolutions, 0 - 1000+
%
%  dim3                     integer             Singleton dimension
%                                               specifier for x,y,z
%                                               components
%
%% Outputs:
%
%  v0                       [N x 3]            Initial Velocity at r0
%                                               and t = 0 to intercept
%                                               rf at t = dt
%
%  vf                       [N x 3]             Final Velocity at rf
%                                               and t = dt
%
%% Revision History:
%  Darin C. Koblick                                              10/21/2025
%  Copyright 2025 Coorbital, Inc.
%% --------------------------- Begin Code Sequence ------------------------
if nargin == 0
        mu = 398600.4418;
        dt = 86400*10;
        r0 = [-13491.4321511262
              -42677.1843763108
              -201.366208254528]';
       rf = [ -4153.6338861922
              12185.7209099644
             -27979.3969414587]';
      dir = -1;
     nmax = 1000;
  [v0,vf] = pumpkyn.pykep.lambert2Body(r0,rf,dt,mu,dir,nmax,2);
  return;
end
%% Flatten inputs:
[r0,fSeq] = pumpkyn.util.fDim(r0,dim3);
       rf = pumpkyn.util.fDim(rf,dim3);
       
%% Call lambert solver:
try
    [v0,vf] = pumpkyn.pykep.pykep_lambert_mex(r0,rf,dt,mu,dir,nmax);
catch
    %Compile mex file, then run:
    mex COPTIMFLAGS="-Ofast -fp:fast" pykep_lambert_mex.cpp ...
        elliptic_orbit.cpp lambert_original.cpp
    [v0,vf] = pumpkyn.pykep.pykep_lambert_mex(r0,rf,dt,mu,dir,nmax);
end
          

%% Reshape outputs:
fSeq.postShift(1) = size(v0,2);
               v0 = pumpkyn.util.eDim(v0,fSeq);
               vf = pumpkyn.util.eDim(vf,fSeq);
end