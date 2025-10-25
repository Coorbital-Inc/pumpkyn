function tmpVar = getConst(varargin)
%% Purpose;
%  Routine which will allow for the user to identify a specific group of
%  constants in which he would like to load into his MATLAB session.
%
%
%% Inputs:
% varargin                     string             Specifies the grouping of
%                                                 constants in which you
%                                                 would like to load
%                                                 If you are not very
%                                                 particular, don't specify
%                                                 this constraint and all
%                                                 constants will load
%
%% Outputs:
%  tmpVar                   struct                Structure containing all
%                                                 of the constants in which
%                                                 the user elected to
%                                                 output at runtime.
%
%% Revision History:
%  Darin C. Koblick                               (c) 2013
%  Copyright 2025 Coorbital, Inc.
%% --------------------- Begin Code Sequence ------------------------------
if nargin == 0
    varargin = {'all'};
end

if any(strcmp(varargin,'astro')) || any(strcmp(varargin,'all'))
%Load the Standard Gravitational Parameters of the Solar System (km^3s^-2)
      tmpVar.Sun.Mu = 132712440018;
    tmpVar.Sun.Mass = 1988500E24;
  tmpVar.Mercury.Mu = 22032;
    tmpVar.Venus.Mu = 324859;
    tmpVar.Earth.Mu = 398600.4418;
tmpVar.EarthMoon.Mu = 4902.8000;     
     tmpVar.Mars.Mu = 42828;        
    tmpVar.Ceres.Mu = 63.1;
  tmpVar.Jupiter.Mu = 126686534;
   tmpVar.Saturn.Mu = 37931187;
   tmpVar.Saturn.Mass = 568.34E24;
   tmpVar.Saturn.Rad = 60268;
   tmpVar.Uranus.Mu = 5793939;
  tmpVar.Neptune.Mu = 6836529;
    tmpVar.Pluto.Mu = 871;
     tmpVar.Eris.Mu = 1108;
%Load the Mean Raidus of the Planets in the Solar System (km) 
    tmpVar.Earth.Rad = 6378.137;         %Vallado
    tmpVar.Earth.Mass =  5.9726e24;      %http://nssdc.gsfc.nasa.gov/planetary/factsheet/earthfact.html
      tmpVar.Sun.Rad = 696000;           %http://nssdc.gsfc.nasa.gov/planetary/factsheet/sunfact.html
     tmpVar.Mars.Rad = 3389.5;           %http://nssdc.gsfc.nasa.gov/planetary/factsheet/marsfact.html
     tmpVar.EarthMoon.Rad = 1738.1;      %http://nssdc.gsfc.nasa.gov/planetary/factsheet/moonfact.html
     tmpVar.EarthMoon.Mass = 0.07346e24; %http://nssdc.gsfc.nasa.gov/planetary/factsheet/moonfact.html
     tmpVar.Venus.Rad = 6051.8;          %https://nssdc.gsfc.nasa.gov/planetary/factsheet/venusfact.html
%Load the Saturn Moon, Titan             %https://nssdc.gsfc.nasa.gov/planetary/factsheet/saturniansatfact.html
    tmpVar.SaturnTitan.Rad = 2575;
   tmpVar.SaturnTitan.Mass = 1345.5*10^20;
     tmpVar.SaturnTitan.Mu =   tmpVar.SaturnTitan.Mass*(6.67384E-11)/1000^3;
    tmpVar.PlutoCharon.Rad = 606;
   tmpVar.PlutoCharon.Mass = 1.586e+21;
     tmpVar.PlutoCharon.Mu = 106.017344589409;
     
     
     
%Load the Blackbody Temperatures in K:
     tmpVar.Sun.T = 5778;                %Effective temperature Solar Atmosphere
                                         %https://radiojove.gsfc.nasa.gov/education/educ/sun/basics/material/sunfacts.htm     
%Load the Rotation Rate of the Earth (rad/sec)     
     tmpVar.Earth.omega = 0.0000729211585530;   %rad/solar sec
         tmpVar.Earth.T = 254.3;        %Blackbody temperature of Earth (K) 
                                        %Source: (NASA Earth Fact Sheet)
    tmpVar.Earth.Albedo = 0.434;        %Visual Geometric Albedo of Earth 
                                        %Source: (NASA Earth Fact Sheet)
                                        %https://nssdc.gsfc.nasa.gov/planetary/factsheet/moonfact.html
    tmpVar.EarthMoon.Albedo = 0.12;     
        tmpVar.Earth.J2 =  0.00108263;  %Source: aerospaceengineering.net      
        tmpVar.Earth.J3 = -0.0000025327;%Source: vallado
        tmpVar.Earth.J4 = -0.0000016196;%Source: vallado
                                        
        tmpVar.AU = 149597871;          %1 AU distance in Kilometers
end

if any(strcmp(varargin,'physics')) || any(strcmp(varargin,'all'))
   tmpVar.Boltz =  5.67051e-8;          %Stefan Boltzman Constant W/(m^2*K^4)
   tmpVar.kBolt =  1.38038426460289e-23;%Boltzmann's Constant (J/K)
   tmpVar.Planck = 6.6260755e-34;       %Planck's Constant (W s^2)
   tmpVar.c      = 299792000.458;       %Speed of Light (m/s)
   tmpVar.uRadPerArcs = 4.84813681;     %Microradians per arcsec
   tmpVar.qElectron = 1.60219E-19;      %Electron charge, C
   tmpVar.deg2ArcSec = 35600;           %3600 arcsec in one deg 
   tmpVar.VmSun =  -26.74;              %Visual Magnitude of the Sun at 1 AU
   tmpVar.VmMoon = -12.74;              %Visual Magnitude of Full Moon at surface of Earth
   tmpVar.VmEarth = -3.86;              %Visual Magnitude of Earth at 1 AU
   tmpVar.G = 6.67384E-11;              %Gravitational Constant m^3 kg^-1 s^-2
   tmpVar.Ratm = 287.058;               %Specific Gas Constant of Dry Air (J/(kg*K))
   tmpVar.R = 8.3144598;                %Univeral Gas Constant (J*mol^-1*K^-1)
   tmpVar.psi2Pa = 6894.75729;          %Convert PSI to Pascals
   tmpVar.g = 9.80665;                  %Earth Gravity Acceleration (Standard) m/s^2 
   tmpVar.in2m = 0.0254;                %Conversion from inches to meters 
   tmpVar.m2in = 1/tmpVar.in2m;         %Conversion from meters to inches
end
end