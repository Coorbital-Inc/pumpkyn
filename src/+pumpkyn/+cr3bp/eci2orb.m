function oev = eci2orb(mu,r,v,dim3)
%% Purpose:
%
% convert eci state vector to six classical orbital
% elements via keplerian orbital elements

%% inputs
%
%  mu             double                        Central body gravitational 
%                                               constant (km**3/sec**2)
%
%  r              [N x 3]                       ECI position vector (km)
%
%  v              [N x 3]                       ECI velocity vector (km/s)
%
%  dim3           integer                       Singleton Dimension
%                                               Specifier
%
%% output
%
% oev             [N x 6]                       Keplerian Orbital Elements
%                                               sma [km]
%                                               ecc
%                                               inc [rad]
%                                               argp [rad]
%                                               raan [rad]
%                                               trua [rad]
%
%% Revision History:
%  Darin C. Koblick                                         (c) 06-10-2023
%% ------------------------- Begin Code Sequence --------------------------
%% Step One: Flatten the input dimensions:
[rv,fSeq] = pumpkyn.util.fDim(cat(dim3,r,v),dim3);
%% Step Two: Call the eci2orbl routine:
      oev = eci2orb1(mu,rv(:,1:3),rv(:,4:6));
%% Step Three: Reshape the outputs:
      oev = pumpkyn.util.eDim(oev,fSeq);
end

function oev = eci2orb1(mu, r, v)

% convert eci state vector to six classical orbital
% elements via keplerian orbital elements

% input

%  mu = central body gravitational constant (km**3/sec**2)
%  r  = eci position vector (kilometers)
%  v  = eci velocity vector (kilometers/second)

% output

%  oev(1) = semimajor axis (kilometers)
%  oev(2) = orbital eccentricity (non-dimensional)
%           (0 <= eccentricity < 1)
%  oev(3) = orbital inclination (radians)
%           (0 <= inclination <= pi)
%  oev(4) = argument of perigee (radians)
%           (0 <= argument of perigee <= 2 pi)
%  oev(5) = right ascension of ascending node (radians)
%           (0 <= raan <= 2 pi)
%  oev(6) = true anomaly (radians)
%           (0 <= true anomaly <= 2 pi)

% Orbital Mechanics with Matlab

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if size(r,1) > 1 && size(r,2) == 3
       oev = NaN(size(r,1),6);
    for tEl=1:size(r,1)
       oev(tEl,:) = eci2orb1(mu,r(tEl,:),v(tEl,:)); 
    end
    return;
end


pi2 = 2.0 * pi;

% position and velocity magnitude

rmag = norm(r);

vmag = norm(v);

% position unit vector

rhat = r / rmag;

% angular momentum vectors

hv = cross(r, v);

hhat = hv / norm(hv);

% eccentricity vector

vtmp = v / mu;

ecc = cross(vtmp, hv);

ecc = ecc - rhat;

% semimajor axis

sma = 1.0 / (2.0 / rmag - vmag * vmag / mu);

p = hhat(1) / (1.0 + hhat(3));

q = -hhat(2) / (1.0 + hhat(3));

const1 = 1.0 / (1.0 + p * p + q * q);

fhat(1) = const1 * (1.0 - p * p + q * q);
fhat(2) = const1 * 2.0 * p * q;
fhat(3) = -const1 * 2.0 * p;

ghat(1) = const1 * 2.0 * p * q;
ghat(2) = const1 * (1.0 + p * p - q * q);
ghat(3) = const1 * 2.0 * q;

h = dot(ecc, ghat);

xk = dot(ecc, fhat);

x1 = dot(r, fhat);

y1 = dot(r, ghat);

% orbital eccentricity

eccm = sqrt(h * h + xk * xk);

% orbital inclination

inc = 2.0 * atan(sqrt(p * p + q * q));

% true longitude

xlambdat = atan3(y1, x1);

% check for equatorial orbit

if (inc > 0.00000001)
    raan = atan3(p, q);
else
    raan = 0.0;
end

% check for circular orbit

if (eccm > 0.00000001)
    argper = mod(atan3(h, xk) - raan, pi2);
else
    argper = 0.0;
end

% true anomaly

tanom = mod(xlambdat - raan - argper, pi2);

% load orbital element vector
oev(1) = sma;
oev(2) = eccm;
oev(3) = inc;
oev(4) = argper;
oev(5) = raan;
oev(6) = tanom;
end

function y = atan3(a,b)
%% Purpose:
%
% four quadrant inverse tangent see:  Orbital Mechanics with MATLAB
%
%% inputs:
%  a = sine of angle
%  b = cosine of angle
%% outputs:
%  y = angle (radians; 0 =< c <= 2 * pi)
%
%% Revision History:
%  Darin C. Koblick                                         (c) 04-29-2018
%% ----------------------- Begin Code Sequence ----------------------------
if nargin == 0
   a = ones(100,1).*pi + randn(100,1);
   b = ones(100,1).*pi + randn(100,1);
   y = atan3(a,b);
   return;
end
a = real(a);
b = real(b);
y = NaN(size(a));
epsilon = 1e-12; %0.0000000001;
pidiv2 = 0.5*pi;
idx = abs(a) < epsilon;
if any(idx)
    y(idx) = (1 - sign(b(idx))).*pidiv2;
end
if any(~idx)
    c = (2 - sign(a(~idx))).*pidiv2;
    at = a(~idx);
    bt = b(~idx);
    eidx = abs(bt) < epsilon;
    y(eidx)  = c(eidx);
    y(~eidx) = c(~eidx) + ...
         sign(at(~eidx)).*sign(bt(~eidx)).* ...
         (abs(atan(at(~eidx)./bt(~eidx))) - pidiv2);
end
end

