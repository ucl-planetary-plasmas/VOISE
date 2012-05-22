function params = getHSTPlanetParams(params)
% function params = getHSTPlanetParams(params)

%
% $Id: getHSTPlanetParams.m,v 1.1 2012/05/22 16:34:59 patrick Exp $
%
% Copyright (c) 2012 Patrick Guio <patrick.guio@gmail.com>
% All Rights Reserved.
%
% This program is free software; you can redistribute it and/or modify it
% under the terms of the GNU General Public License as published by the
% Free Software Foundation; either version 3 of the License, or (at your
% option) any later version.
%
% This program is distributed in the hope that it will be useful, but
% WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General
% Public License for more details.
%
% You should have received a copy of the GNU General Public License
% along with this program. If not, see <http://www.gnu.org/licenses/>.

HST = params.HST;

isJup = strfind(HST.TARGNAME, 'JUP');
isSat = strfind(HST.TARGNAME, 'SAT');

if ~isempty(isJup) & isJup,
  planet = 'jupiter';
elseif ~isempty(isSat) & isSat,
  planet = 'saturn';
else
  planet = lower(HST.TARGNAME);
end


%PLANET = upper(planet);
%Planet = [PLANET(1) planet(2:end)];
IAU_PLANET = ['IAU_' upper(planet)];

degPerRad = cspice_dpr;
radPerDeg = cspice_rpd;

radPerArcsec = radPerDeg/3600;
arcsecPerRad = degPerRad*3600;

% generic kernel path
spiceKernelsPath = getSpiceGenericKernelsPath();

% Load a leapseconds kernel.
% ftp://naif.jpl.nasa.gov/pub/naif/generic_kernels/lsk/
cspice_furnsh([spiceKernelsPath 'naif0010.tls']);

% Load planetary ephemeris 
% ftp://naif.jpl.nasa.gov/pub/naif/generic_kernels/spk/planets/
cspice_furnsh([spiceKernelsPath 'de421.bsp']);

% Load satellite ephemeris 
% ftp://naif.jpl.nasa.gov/pub/naif/generic_kernels/spk/satellites/
switch planet,
  case 'jupiter'
    cspice_furnsh([spiceKernelsPath 'jup291.bsp']);
  case 'saturn'
    cspice_furnsh([spiceKernelsPath 'sat353.bsp']);
  case 'uranus'
    cspice_furnsh([spiceKernelsPath 'ura095.bsp']);
end

% Load orientation data for planets, natural 
% satellites, the Sun, and selected asteroids
% ftp://naif.jpl.nasa.gov/pub/naif/generic_kernels/pck/
cspice_furnsh([spiceKernelsPath 'pck00010.tpc']);


% Convert string date to cspice
et = cspice_str2et([HST.TDATEOBS ' ' HST.TTIMEOBS]);
% Centre time half way exposure time
%et = et + HST.EXPTIME/2;

% Get position of planet with respect to Earth
target   = planet;
frame    = 'IAU_EARTH';
abcorr   = 'NONE'; 
abcorr   = 'LT';   % Correct for one-way light time
%abcorr   = 'LT+S'; % Correct for one-way light time and stellar aberration
observer = 'EARTH';
% Look up the 'state' vectors and light time values 'ltime'  
% corresponding to the vector of input ephemeris time 'et'.
[state , ltime] = cspice_spkezr(target, et, frame, abcorr, observer);
planetposn = state(1:3)
planetdist = cspice_vnorm(planetposn);

% planet ra and dec in J2000 ref frame
[state , ltime] = cspice_spkezr(target, et, 'J2000', abcorr, observer);
[planetrange,planetra,planetdec] = cspice_recrad(state(1:3));
params.planet.ra = planetra*degPerRad;
params.planet.rec = planetdec*degPerRad;

% Convert the RA/DEC values to radians.
if 0
ra  = HST.RA_TARG * radPerDeg;
dec = HST.DEC_TARG * radPerDeg;
elseif 0
ra  = HST.CRVAL1 * radPerDeg;
dec = HST.CRVAL2 * radPerDeg;
else
ra  = HST.RA_APER * radPerDeg;
dec = HST.DEC_APER * radPerDeg;
end
fprintf(1,'planet    ra, dec = %.8f, %.8f\n',params.planet.ra,params.planet.rec);
fprintf(1,'RA_TARG, DEC_TARG = %.8f, %.8f\n', HST.RA_TARG, HST.DEC_TARG);
fprintf(1,' CRVAL1,   CRVAL2 = %.8f, %.8f\n', HST.CRVAL1, HST.CRVAL2);
fprintf(1,'RA_APER, DEC_APER = %.8f, %.8f\n', HST.RA_TARG, HST.DEC_TARG);

pause

% Convert the angular description of the unit vectors to rectangular
% coordinates.
targetposn = cspice_radrec(planetdist, ra, dec);


% Retrieve the transformation matrix from frames J2000 to IAU_EARTH.
J200toEarth = cspice_pxform('J2000', 'IAU_EARTH', et);

% transform rectangular coordinates into Earth frame
targetposn = J200toEarth*targetposn;

fprintf(1,'planetposn %12.6g, %12.6g %12.6g\n', planetposn);
fprintf(1,'targetposn %12.6g, %12.6g %12.6g\n', targetposn);

target2planetangle = 3600*acosd(dot(cspice_vhat(planetposn),cspice_vhat(targetposn)))

% normalised radial vector from Earth toward target (line of sight) 
z = cspice_vhat(targetposn);
% celestial north is y in Earth frame of reference
north = [0;0;1];
% y is projection of celestial north onto the image plane
% perpendicular to the line of sight
%cspice_vnorm(north-dot(z,north)*z)
y = cspice_vhat(north-dot(z,north)*z);
% angle between celestial north and its projection onto the image plane
tprojnorth = acosd(dot(north,cspice_vhat(north-dot(z,north)*z)))
zdoty = dot(z,y), % scalar product should be zero
% x such that (x,y,z) is direct
x = cross(y,z);

% target position (xt,yt) should be zero
xt = dot(targetposn,x);
yt = dot(targetposn,y);
% km to pixel
xt = atan2(xt,planetdist)*arcsecPerRad/HST.PLATESC
yt = atan2(yt,planetdist)*arcsecPerRad/HST.PLATESC

%pause

% planet centre (xpc,ypc) in the image plane in km
% with respect to (0,0) that corresponds to targetposn
xpc = dot(planetposn,x);
ypc = dot(planetposn,y);
% km to pixel
xpc = atan2(xpc,planetdist)*arcsecPerRad/HST.PLATESC
ypc = atan2(ypc,planetdist)*arcsecPerRad/HST.PLATESC

% centre of image
if 0
ximc = size(HST.W,2)/2;
yimc = size(HST.W,1)/2;
else
ximc = HST.CRPIX1;
yimc = HST.CRPIX2;
end

% rotation to position angle of image y axis (deg. e of n)
% i.e. rotated by ORIENTAT corresponding to the position
% angle of image y axis (in degrees East of North)
[xpcrot,ypcrot] = rot(HST.ORIENTAT,xpc,ypc);

%pause

% Retrieve the transformation matrix from frame of the planet to IAU_EARTH.
planet2Earth = cspice_pxform(IAU_PLANET, 'IAU_EARTH', et);
% Jupiter axis in Earth frame
planetaxis = planet2Earth*[0;0;1];
% angle between rotation axis and projection onto the image plane
tprojplanetaxis = acosd(dot(planetaxis,cspice_vhat(planetaxis-dot(z,planetaxis)*z)))
% projection onto the plane of the image
planetaxis = cspice_vhat(planetaxis-dot(z,planetaxis)*z)
xplanetaxis = dot(planetaxis,x);
yplanetaxis = dot(planetaxis,y);
[xplanetaxis,yplanetaxis] = rot(HST.ORIENTAT,xplanetaxis,yplanetaxis);
planetaxis = [xplanetaxis,yplanetaxis];

% get planet radii
radii = cspice_bodvrd(planet,'RADII',3);
a = radii(1);
b = radii(3);
e = sqrt(a^2-b^2)/a;

if 0
% correction for projection on the plane
a = a*cosd(tprojplanetaxis);
b = b*cosd(tprojplanetaxis);
end

% projection of semi-major and semi minor axis in pixels
a = atan2(a,planetdist)*arcsecPerRad/HST.PLATESC
b = atan2(b,planetdist)*arcsecPerRad/HST.PLATESC
% tilt in the plane of the image with respect to x axis
tilt = atan2(planetaxis(2),planetaxis(1))*degPerRad

% plot an ellispe with planet parameters
theta = linspace(0,2*pi,100);
ellx = a*cos(theta);
elly = b*sin(theta);

xAxis = [-1;1]*planetaxis(1)*b;
yAxis = [-1;1]*planetaxis(2)*b;

% rotate ellipse of tilt 
[ellx,elly] = rot(tilt,ellx,elly);

%[xAxis,yAxis] = rot(tilt,xAxis,yAxis);

W = params.W;
plot(ellx+xpc+ximc,elly+ypc+yimc,'k-',xAxis+xpc+ximc,yAxis+ypc+yimc,'k-')
hold on
W(params.W<=.2) = .2;
X = [1:size(W,2)];
Y = [1:size(W,1)];
imagesc(X,Y,log10(W))
axis xy
axis equal
plot(ximc,yimc,'+','markersize',10)
plot([ximc,xpc+ximc],[yimc,ypc+yimc],'-x','markersize',10)
plot(ellx+xpc+ximc,elly+ypc+yimc,'k-',xAxis+xpc+ximc,yAxis+ypc+yimc,'k-')
hold off
%pause


% rotate of ORIENTAT
[ellx,elly] = rot(HST.ORIENTAT,ellx,elly);

[xAxis,yAxis] = rot(HST.ORIENTAT,xAxis,yAxis);

hold on
plot([ximc,xpcrot+ximc],[yimc,ypcrot+yimc],'-x','markersize',10)
plot(ellx+xpcrot+ximc,elly+ypcrot+yimc,'-',xAxis+xpcrot+ximc,yAxis+ypcrot+yimc)
hold off


%  It's always good form to unload kernels after use,
%  particularly in MATLAB due to data persistence.
cspice_kclear


function [X,Y] = rot(alpha,x,y)

cosa = cosd(alpha);
sina = sind(alpha);

X = x*cosa - y*sina;
Y = x*sina + y*cosa;

