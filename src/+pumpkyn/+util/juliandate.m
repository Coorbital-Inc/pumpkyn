function jd = juliandate( varargin ) 
%% Purpose:
% This routine is a replacement for the MATLAB juliandate routine which
% will compute the julian date given strings/datevectors in the same format
% that datenum would compute.
%
% A Julian date is the number of days and fractional days since noon on 
% November 24, 4713 BCE in the proleptic Gregorian calendar, 
% or January 1, 4713 BCE in the proleptic Julian calendar.
%
%% Examples:
%		d = '12/24/1984';
%		t = 725000.00;
%		c = datevec(d) or c = datevec(t) produce c = [1984 12 24 0 0 0].
%		[y,m,d,h,mi,s] = datevec(d) returns y=1984, m=12, d=24, h=0, mi=0, s=0.
%		c = datevec('5/6/03') produces c = [2003 5 6 0 0 0] until 2054.
%		c = datevec('5/6/03',1900) produces c = [1903 5 6 0 0 0].
%		c = datevec('19.05.2000','dd.mm.yyyy') produces c = [2000 5 19 0 0 0].
%
%% Outputs:
%  jd
%% Revision History:
%  Darin C. Koblick                                         (c) 03-11-2016
%% ---------------------- Begin Code Sequence -----------------------------

[year,month,day,hour,min,sec] = datevec(datenum(varargin{:}));

for k = length(month):-1:1
    if ( month(k) <= 2 ) % January & February
        year(k)  = year(k) - 1.0;
        month(k) = month(k) + 12.0;
    end
end

jd = floor( 365.25*(year + 4716.0)) + floor( 30.6001*( month + 1.0)) + 2.0 - ...
    floor( year/100.0 ) + floor( floor( year/100.0 )/4.0 ) + day - 1524.5 + ...
    (hour + min/60 + sec/3600)/24;
