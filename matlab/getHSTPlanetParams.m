function params = getHSTPlanetParams(params)
% function params = getHSTPlanetParams(params)

%
% $Id: getHSTPlanetParams.m,v 1.6 2012/06/13 16:58:47 patrick Exp $
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

% load necessary SPICE kernels
planet = setUpSpice4Planet(HST);

IAU_PLANET = ['IAU_' upper(planet.name)];

% conversion factors
degPerRad = cspice_dpr;
radPerDeg = cspice_rpd;
radPerArcsec = radPerDeg/3600;
arcsecPerRad = degPerRad*3600;

radPerPixel = arcsecPerRad/HST.PLATESC;

% Convert string date to double precision 
et = cspice_str2et([HST.TDATEOBS ' ' HST.TTIMEOBS]);

% flag for calculation in J2000
if 0,
isJ2000 = true;
else
isJ2000 = false;
end

% Get position of planet with respect to Earth
target   = planet.name;
if isJ2000,
frame    = 'J2000'; % Earth inertial frame
else
frame    = 'IAU_EARTH'; % Earth fixed frame
end
abcorr   = 'LT';   % Correct for one-way light time
observer = 'EARTH';
% Look up the 'state' vectors and light time values 'lt'  
% corresponding to the vector of input ephemeris time 'et'.
[state , lt] = cspice_spkezr(target, et, frame, abcorr, observer);
planetposn = state(1:3);
planetdist = cspice_vnorm(planetposn);

% Get planet ra/dec requires state in inertial J2000
[state , lt] = cspice_spkezr(target, et, 'J2000', 'LT', observer);
[range,ra,dec] = cspice_recrad(state(1:3,:));
planet.r = range;
% rad to deg
planet.ra = ra*degPerRad; 
planet.dec = dec*degPerRad;

fprintf(1,'planet    ra, dec = %12.6f,%12.6f deg\n',[planet.ra(:)';planet.dec(:)']);

% reference pixel image coordinates 
rpx = HST.CRPIX1;
rpy = HST.CRPIX2;
% reference pixel ra/dec coordinates (deg)
rpra = HST.CRVAL1;
rpdec = HST.CRVAL2;
% inverse matrix to transform from world to pixel coordinates
iCD = HST.iCD;
% planet world coordinates to pixel coordinates
[xpc,ypc] = getHSTradec2pixel(HST,planet.ra,planet.dec);

% Convert reference pixel ra/dec (target direction) to rect coordinates
refpixposn = cspice_radrec(planetdist, rpra*radPerDeg, rpdec*radPerDeg);
% Retrieve the transformation matrix from J2000 frame to IAU_EARTH frame
J2000toEarth = cspice_pxform('J2000', 'IAU_EARTH', et);
EarthtoJ2000 = cspice_pxform('IAU_EARTH', 'J2000', et);
% transform rectangular inertial coordinates into Earth frame
if ~isJ2000,
refpixposn = J2000toEarth*refpixposn;
end

fprintf(1,'planetposn        = %12.6g, %12.6g, %12.6g km\n', planetposn);
fprintf(1,'refpixposn        = %12.6g, %12.6g, %12.6g km\n', refpixposn);

refpix2planetangle = acosd(dot(cspice_vhat(planetposn),cspice_vhat(refpixposn)));
fprintf(1,'angle(planet,refpix)           = %f deg\n', refpix2planetangle);
fprintf(1,'|planet(ra,dec)-refpix(ra,dec) = %f deg\n', ...
norm([planet.ra-rpra;planet.dec-rpdec],2));

% zhat normalised vector from ref pixel pointing toward Earth in Earth ref
zhat = -cspice_vhat(refpixposn);
% get celestial north in Earth ref
if isJ2000,
north = EarthtoJ2000*[0;0;1];
else
north = [0;0;1];
end
% yhat is the projection of celestial north onto the image plane
% (the plane perpendicular to the line of sight defined by zhat)
% rotated by ORIENTAT - position angle of image y axis (deg. e of n)
yhat = cspice_vhat(north-dot(zhat,north)*zhat);
% rotate by ORIENTAT corresponding to the position angle of
% image y axis (deg. e of n)
yhatr = rot3d(zhat,HST.ORIENTAT,yhat);
%yhat = yhatr;
% angle between celestial north and its projection onto the image plane
%tprojnorth = acosd(dot(north,cspice_vhat(north-dot(zhat,north)*zhat)))
%zdoty = dot(zhat,yhat), % scalar product should be zero
% xhat such that (xhat,yhat,zhat) is direct
xhat = cross(yhat,zhat);
xhatr = cross(yhatr,zhat);
if 0
fprintf(1,'angle(y/xhat,y/xhatr)=%.2f,%.2f\n', ...
        [acosd(dot(yhat,yhatr)), acosd(dot(xhat,xhatr))]);
fprintf(1,'xhat(r).yhat(r)      =%.2g,%.2g\n',[dot(xhat,yhat), dot(xhatr,yhatr)])
end
% rotated axes
xhat = xhatr; yhat = yhatr;

% ref pixel position (xrp,yrp) should be zero
xrp = atan2(dot(refpixposn,xhat), planetdist)*radPerPixel;
yrp = atan2(dot(refpixposn,yhat), planetdist)*radPerPixel;
fprintf(1,'xrp, yrp          = %12.6g, %12.6g pixel\n', xrp, yrp);

if 1
% planet centre (Xpc,Ypc) in the image plane in km
% with respect to (0,0) that corresponds to refpixposn
Xpc = atan2(dot(planetposn,xhat), planetdist)*radPerPixel + rpx;
Ypc = atan2(dot(planetposn,yhat), planetdist)*radPerPixel + rpy; 
end

pc = [xpc-HST.CRPIX1, ypc-HST.CRPIX2];
PC = [Xpc-HST.CRPIX1, Ypc-HST.CRPIX2];
fprintf(1,'xpc, ypc, |pc|    = %12.6f, %12.6f, %12.6f pixel\n', pc, norm(pc));
fprintf(1,'xPC, yPC, |PC|    = %12.6f, %12.6f, %12.6f pixel\n', PC, norm(PC));
fprintf(1,'angle(pc,PC)      = %12.2f deg\n', acosd(dot(pc/norm(pc),PC/norm(PC))));


% get planet radii
radii = cspice_bodvrd(planet.name,'RADII',3);
a = radii(1); % equatorial radius (1 and 2)
b = radii(3); % polar radius
e = sqrt(a^2-b^2)/a; % excentricity

% if Saturn get rings parameters
if strcmp(planet.name,'saturn'),
  % A Ring, Encke Gap, Cassini Division, B Ring
  rings = {'RING1','RING1_1','RING2','RING3'};
	ringSpecs = zeros(5,length(rings));
	% Ring geometry is defined in the form of one set of R1, R2, Z1, Z2, OD where
	% R1 and R2 are inner and outer radii of the ring (in km)
	% Z1 and Z2 are the vertical heights of the ring at R1 and R2 (also in km, 
	% equal to one-half of the total thickness of the ring) 
	% OD is the average optical depth of the ring sub-segment/gap across R1 to R2.
	for i=1:length(rings),
   ringSpecs(:,i) = cspice_bodvrd(planet.name,rings{i},5);
	end
end

% Retrieve the transformation matrix from frame of the planet to IAU_EARTH.
if 1, 
% one-way light time corrected epoch
planet2Earth = cspice_pxform(IAU_PLANET, 'IAU_EARTH', et-lt);
else
planet2Earth = cspice_pxform(IAU_PLANET, 'IAU_EARTH', et);
end
% planet axis direction (pointing north) in Earth frame (or inertial J2000)
if ~isJ2000,
planetaxis = planet2Earth*[0;0;1];
else
planetaxis = EarthtoJ2000*(planet2Earth*[0;0;1]);
end


% angle between planet axis and projection onto the image plane
tprojplanetaxis = acosd(dot(planetaxis,cspice_vhat(planetaxis-dot(zhat,planetaxis)*zhat)));
fprintf(1,'tprojplanetaxis      = %12.2g deg\n', tprojplanetaxis);
fprintf(1,'cos(tprojplanetaxis) = %12.2g\n', cosd(tprojplanetaxis)),
% projection of the planet axis onto the plane of the image
planetaxis = cspice_vhat(planetaxis-dot(zhat,planetaxis)*zhat);
fprintf(1,'planetaxis (Earth r) = %12.6f, %12.6f, %12.6f\n', planetaxis);
xplanetaxis = dot(planetaxis,xhat);
yplanetaxis = dot(planetaxis,yhat);
planetaxis3d = planetaxis;
planetaxis = [xplanetaxis,yplanetaxis];
fprintf(1,'planetaxis (image)   = %12.6f, %12.6f\n', planetaxis);

spheroidGrid(planetaxis3d,a,b,xhat,yhat,zhat);


% poles computed in world (ra/dec) and transformed to pixel coordinates
if ~isJ2000,
planetnorth = EarthtoJ2000*(planetposn+planet2Earth*[0;0;b]);
planetsouth = EarthtoJ2000*(planetposn+planet2Earth*[0;0;-b]);
else
planetnorth = planetposn+EarthtoJ2000*(planet2Earth*[0;0;b]);
planetsouth = planetposn+EarthtoJ2000*(planet2Earth*[0;0;-b]);
end
[northrange,northra,northdec] = cspice_recrad(planetnorth);
[southrange,southra,southdec] = cspice_recrad(planetsouth);
[xn,yn] = getHSTradec2pixel(HST,northra*degPerRad,northdec*degPerRad);
[xs,ys] = getHSTradec2pixel(HST,southra*degPerRad,southdec*degPerRad);

% equator computed in world (ra/dec) and transformed to pixel coordinates,
lat = 0;
lon = linspace(0,2*pi,100);
planeteq = zeros(3,length(lon));
for i=1:length(lon);
  posn = [a*cos(lat)*cos(lon(i));a*cos(lat)*sin(lon(i));b*sin(lat)];
	if ~isJ2000,
	planeteq(:,i) = EarthtoJ2000*(planetposn+planet2Earth*posn);
	else
	planeteq(:,i) = planetposn+EarthtoJ2000*(planet2Earth*posn);
	end
end
[eqrange,eqra,eqdec] = cspice_recrad(planeteq);
[xeq,yeq] = getHSTradec2pixel(HST,eqra*degPerRad,eqdec*degPerRad);

% meridian computed in world (ra/dec) and transformed to pixel coordinates
lat = linspace(-pi/2,pi/2,100);
lon = 0;
planetmerid = zeros(3,length(lat));
for i=1:length(lat);
  posn = [a*cos(lat(i))*cos(lon(1));a*cos(lat(1))*sin(lon(1));b*sin(lat(i))];
	if ~isJ2000,
  planetmerid(:,i) = EarthtoJ2000*(planetposn+planet2Earth*posn);
	else
  planetmerid(:,i) = planetposn+EarthtoJ2000*(planet2Earth*posn);
	end
end
[meridrange,meridra,meriddec] = cspice_recrad(planetmerid);
[xmerid,ymerid] = getHSTradec2pixel(HST,meridra*degPerRad,meriddec*degPerRad);


if strcmp(planet.name,'saturn'),
% Saturn's ring computed in world (ra/dec) and transformed to pixel coordinates,
lon = linspace(0,2*pi,100);
for j = 1:length(rings),
rmn{j} =  zeros(3,length(lon));
rmx{j} =  zeros(3,length(lon));
for i=1:length(lon);
  posn = [ringSpecs(1,j)*cos(lon(i));ringSpecs(1,j)*sin(lon(i));0];
	if ~isJ2000,
	rmn{j}(:,i) = EarthtoJ2000*(planetposn+planet2Earth*posn);
	else
	rmn{j}(:,i) = planetposn+EarthtoJ2000*(planet2Earth*posn);
	end
  posn = [ringSpecs(2,j)*cos(lon(i));ringSpecs(2,j)*sin(lon(i));0];
	if ~isJ2000,
	rmx{j}(:,i) = EarthtoJ2000*(planetposn+planet2Earth*posn);
	else
	rmx{j}(:,i) = planetposn+EarthtoJ2000*(planet2Earth*posn);
	end
end
[rmnr{j},rmnra{j},rmndec{j}] = cspice_recrad(rmn{j});
[xrmn{j},yrmn{j}] = getHSTradec2pixel(HST,rmnra{j}*degPerRad,rmndec{j}*degPerRad);
[rmxr{j},rmxra{j},rmxdec{j}] = cspice_recrad(rmx{j});
[xrmx{j},yrmx{j}] = getHSTradec2pixel(HST,rmxra{j}*degPerRad,rmxdec{j}*degPerRad);
end
end

if 1,
% correction for projection on the plane
a = a*cosd(tprojplanetaxis);
b = b*cosd(tprojplanetaxis);
end

% projection of semi-major and semi minor axis in pixels
a = atan2(a,planetdist)*radPerPixel;
b = atan2(b,planetdist)*radPerPixel;
% tilt in the plane of the image with respect to x axis
tilt = atan2(planetaxis(2),planetaxis(1))*degPerRad;

fprintf(1,'planetary disc  a, b = %12.1f, %12.1f pixel, tilt = %5.1f deg\n',...
a,b,tilt);

if strcmp(planet.name,'saturn'),
% correction for projection on the plane
%ARmin = ARmin*cosd(tprojplanetaxis);
end

% plot an ellipse with planet parameters, i.e. a, b, tilt
theta = linspace(0,360,100);
% ellipse with axis along x axis
ellx = b*sind(theta);
elly = a*cosd(theta);
% rotation by tilt
ell = rot2d(tilt,[ellx;elly]);
ellx = ell(1,:);
elly = ell(2,:);

opts = {'fontsize',12,'fontweight','normal'}; %,'color','black'};

% north-south pole axis
xAxis = [1;-1]*planetaxis(1)*b;
yAxis = [1;-1]*planetaxis(2)*b;

if 0
plot(ellx,elly,xAxis,yAxis)
text(xAxis(1),yAxis(1),'NAxis',opts{:});
text(xAxis(2),yAxis(2),'SAxis',opts{:});
axis equal
pause
end

figure
% load color values
h=plot(ones(10));
for i=1:length(h),
  colors{i} = get(h(i),'color');
end
clear h
close


W = params.W;
clf
X = [1:size(W,2)];
Y = [1:size(W,1)];
imagesc(X,Y,log10(abs(W)))
axis xy
hold on
axis auto
plot(rpx,rpy,'b+','markersize',10), text(rpx,rpy,'O',opts{:})

% using HST transformation from/to world to/from pixel coordinates and SPICE
plot(xpc,ypc,'bx','markersize',10), text(xpc,ypc,'(xpc,ypc)',opts{:})
plot([xn;xs],[yn;ys],'bo-','markersize',10), 
text([xn;xs],[yn;ys],['N';'S'],opts{:},'color','b'),
plot(xeq,yeq,'bo-',xmerid,ymerid,'bo-');
if strcmp(planet.name,'saturn'),
for i=1:length(rings),
plot(xrmn{i},yrmn{i},'o-','color',colors{i});
plot(xrmx{i},yrmx{i},'o-','color',colors{i});
end
end

% using HST ref pixel, y-axis orientation, plate scale and SPICE
plot(Xpc,Ypc,'kx','markersize',10), text(Xpc,Ypc,'(Xpc,Ypc)',opts{:})
plot(ellx+Xpc,elly+Ypc,'k-',xAxis+Xpc,yAxis+Ypc,'g-')
text(xAxis+Xpc,yAxis+Ypc,['NAxis';'SAxis'],opts{:},'color','k')

hold off
axis equal
%pause

% Embed planet structure into params
params.planet = planet;

%  It's always good form to unload kernels after use,
%  particularly in MATLAB due to data persistence.
cspice_kclear


function P = rot2d(alpha,P)

% P is a 2xN vector

x = P(1,:);
y = P(2,:);

ca = cosd(alpha);
sa = sind(alpha);

X = x*ca - y*sa;
Y = x*sa + y*ca;

P(1,:) = X;
P(2,:) = Y;


function P = rot3d(u,alpha,P)

% from http://en.wikipedia.org/wiki/Rotation_matrix

% P is a 3xN vector
x = P(1,:);
y = P(2,:);
z = P(3,:);

ca = cosd(alpha);
sa = sind(alpha);

if 0,
uxx = u(1)*u(1);
uyy = u(2)*u(2);
uzz = u(3)*u(3);
uxy = u(1)*u(2);
uxz = u(1)*u(3);
uyz = u(2)*u(3);;

X = (ca+uxx*(1-ca))*x      + (uxy*(1-ca)-u(3)*sa)*y + (uxz*(1-ca)+u(2)*sa)*z;
Y = (uxy*(1-ca)+u(3)*sa)*x + (ca+uyy*(1-ca))*y      + (uyz*(1-ca)-u(1)*sa)*z;
Z = (uxz*(1-ca)-u(2)*sa)*x + (uyz*(1-ca)+u(1)*sa)*y + (ca+uzz*(1-ca))*z;

else

% identity matrix
I = eye(3,3);
% cross product matrix
u_x = [0, -u(3), u(2); u(3), 0, -u(1); -u(2), u(1), 0];
% tensor product 
uxu = u(:)*u(:)';
% Rotation matrix
R = I*ca + sa*u_x+(1-ca)*uxu;

X = R(1,1)*x + R(1,2)*y + R(1,3)*z;
Y = R(2,1)*x + R(2,2)*y + R(2,3)*z;
Z = R(3,1)*x + R(3,2)*y + R(3,3)*z;

end

P(1,:) = X;
P(2,:) = Y;
P(3,:) = Z;
