function plotPlanetGrid(planet,params,pc,epoch,CML,psi,orientat,PIXSIZE)
% function plotPlanetGrid(planet,params,pc,epoch,CML,psi,orientat,PIXSIZE)

%
% $Id: plotPlanetGrid.m,v 1.3 2015/09/26 19:28:29 patrick Exp $
%
% Copyright (c) 2009 
% Patrick Guio <p.guio@ucl.ac.uk>
%
% All Rights Reserved.
%
% This program is free software; you can redistribute it and/or modify it
% under the terms of the GNU General Public License as published by the
% Free Software Foundation; either version 2.  of the License, or (at your
% option) any later version.
%
% This program is distributed in the hope that it will be useful, but
% WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General
% Public License for more details.
%

%subplot(224)
subplot(111)

% Latitude and longitude spacing for plots [deg]
dlat = 10/2; 
dlon = 20/2;

% Pixel coordinates for planet centre
PCX = pc(1); 
PCY = pc(2);
fprintf(1,'Planet center %8.2f, %8.2f [pixels]\n', PCX,PCY);

% calculation of the subEarth and subSolar latitudes and longitudes
[ss,se] = computePlanetAxis(planet,epoch);
se.CML = CML;
se.psi = psi;

% conversion factor
d2r = cspice_rpd;
r2d = cspice_dpr;

% km to pixel
km2pix = (1/se.distkm)*r2d*3600/PIXSIZE;
% pixel to km -- rad to arcsec is (180/pi)*3600
pix2km = PIXSIZE/3600*d2r*se.distkm;

if length(pc)==2, % no information about ellipsoid
  [a,b,e,f] = getPlanetGeometry(planet);
	% limb eccentricity 
  eL = e*cos(pi/180*se.lat);
	% projected semi minor axis
  bp = a * sqrt(1-eL^2);
	semiMaj_km = a;
	ecc = e;
else, % information about ellispoide provide in pixel
  a = pc(3)*pix2km;
  if length(pc)==3,
    bp = a;
  else,
    bp = pc(4)*pix2km;
  end
  semiMaj_km = a;
  eL = sqrt(1-(bp/a)^2);
	ecc = eL/cos(pi/180*se.lat);
end
fprintf(1,'Semi-major/min axes = %9.4f pixels/%9.4f pixels\n',[semiMaj_km,bp]*km2pix);
fprintf(1,'Semi-major axis/Ecc = %9.0f km    /%9.5f\n',semiMaj_km,ecc);


if ~isempty(params),
% Draw the image data, on a scale of arcsec
% Note that PIXSIZE is the size of one side of a square pixel in arcsec
x = PIXSIZE*(params.x-PCX);
y = PIXSIZE*(params.y-PCY);

imagesc(x, y, params.Wo);
axis xy;
axis equal
axis tight
xlabel('x [arcsec]')
ylabel('y [arcsec]')
%title(sprintf('[%.2f,%.2f]x[%.2f,%.2f] [arcsec]',min(x),max(x),min(y),max(y)));
%title(sprintf('%.2fx%.2f [arcsec]',max(x)-min(x),max(y)-min(y)));
title(sprintf('%s %s',[upper(planet(1)),lower(planet(2:end))],...
      datestr(datenum(epoch,'yyyy mm dd HH MM SS'))));
tic
drawPlanetGrid(planet,epoch,ss,se,orientat,dlat,dlon,semiMaj_km,ecc);
toc
else
tic
drawPlanetGrid(planet,epoch,ss,se,orientat,dlat,dlon,semiMaj_km,ecc);
toc
axis xy;
axis equal
axis tight
xlabel('x [arcsec]')
ylabel('y [arcsec]')
title(sprintf('%s %s',[upper(planet(1)),lower(planet(2:end))],...
      datestr(datenum(epoch,'yyyy mm dd HH MM SS'))));
end



function drawPlanetGrid(planet,epoch,ss,se,orientat,dlat,dlon,semiMaj_km,ecc)

% Get difference between (right-handed) sub-Earth and sub-solar longitudes.
% Although these two quantities will vary on time scale of the planet's
% rotation, their *difference* varies much more slowly. 
% Note that a *positive* value means that the Sun central meridian is
% situated ahead of the Earth's (in the direction of planetary rotation)
ssedlon = ss.lon - se.lon;

% Sub-Earth colatitude in radians [0,pi] -> (NP,SP)
theobs = pi/2 - se.lat*pi/180;
% CML in a right-handed system [0,2pi]
phiobs = 2*pi - se.CML*pi/180;

% Sub-solar colatitude in radians
thesun = pi/2 - ss.lat*pi/180;

% Scaling factor to convert from km to arcsec on the plane of 
% the Earth observer's sky
km2asec = (1/se.distkm)*(180/pi)*3600;

% Calculate the viewing angle
lineOfSight = [sin(theobs)*cos(phiobs);sin(theobs)*sin(phiobs);cos(theobs)];

% http://www.imcce.fr/en/ephemerides/formulaire/form_ephephys.php
% gives 254.31 deg 
alpha = se.psi-orientat;
fprintf(1,'psi %f orientat %f alpha %f (deg)\n',se.psi,orientat,alpha);

hold on

if 0,
A = viewmtx(90-phiobs*180/pi,90-theobs*180/pi);
A(4,4) = 1/km2asec;
[X,Y,Z]=ellipsoid(0,0,0,semiMaj_km,semiMaj_km,semiMaj_km*sqrt(1-ecc^2));
[m,n] = size(X);
x4d = [X(:),Y(:),Z(:),ones(m*n,1)]';
x3d = A*x4d;
x2 = zeros(m,n); y2 = zeros(m,n); z2 = zeros(m,n);
[min(x3d(4,:)),max(x3d(4,:))]
x2(:) = x3d(1,:)./x3d(4,:);
y2(:) = x3d(2,:)./x3d(4,:);
z2(:) = x3d(3,:)./x3d(4,:);
mesh(x2,y2,z2,zeros(size(x2)));hidden off, alpha(1), view(0,90)
pause
f=figure;
subplot(211), mesh(x2,y2,z2,zeros(size(x2))), view(0,90)
subplot(212), plot(x2(z2>=0),y2(z2>=0))
pause
close(f);
end

% Grid curves of constant latitude
for the = ([dlat:dlat:180-dlat])*pi/180,
%for the = ([90+se.lat])*pi/180,
	phi = linspace(0,2*pi,fix(150*sin(the))); 
	[r,x,y,z,xsky,ysky,zsky] = spherical2Sky(semiMaj_km,ecc, ...
	                                    the,phi,theobs,phiobs,km2asec);
  if 0,
  cos_vang = lineOfSight(1)*x+lineOfSight(2)*y+lineOfSight(3)*z;
	flag_vang = (cos_vang > 0);
	else
	flag_vang = (zsky > 0);
	end
	Xv = xsky; Xv(~flag_vang) = NaN;
	Yv = ysky; Yv(~flag_vang) = NaN;
	% rotation of sky plane
	[Xv,Yv] = Rotate(alpha,Xv,Yv);
  plot(Xv,Yv,'k-','LineWidth',.5);
	if 0, % hidden lines
	Xh = xsky; Xh(flag_vang) = NaN;
	Yh = ysky; Yh(flag_vang) = NaN;
  [Xh,Yh] = Rotate(alpha,Xh,Yh);
  plot(Xh,Yh,'k--','LineWidth',.5);
	end
end

% Grid curves of constant longitude
the = linspace(0,pi,50);
for phi = [0:dlon:360-dlon]*pi/180,

	[r,x,y,z,xsky,ysky,zsky] = spherical2Sky(semiMaj_km,ecc, ...
	                                    the,phi,theobs,phiobs,km2asec);
  if 0,
  cos_vang = lineOfSight(1)*x+lineOfSight(2)*y+lineOfSight(3)*z;
  flag_vang = (cos_vang > 0);
  else
	flag_vang = (zsky > 0);
	end
	Xv = xsky; Xv(~flag_vang) = NaN;
	Yv = ysky; Yv(~flag_vang) = NaN;
	% rotation of sky plane
	[Xv,Yv] = Rotate(alpha,Xv,Yv);
	if phi==0, % CML
  plot(Xv,Yv,'b-','LineWidth',2);
	else,
	plot(Xv,Yv,'k-','LineWidth',1.);
	end
	if 0, % hidden lines
	Xh = xsky; Xh(flag_vang) = NaN;
	Yh = ysky; Yh(flag_vang) = NaN;
  [Xh,Yh] = Rotate(alpha,Xh,Yh);
	if phi==0, % CML
  plot(Xv,Yv,'b--','LineWidth',2);
	else
  plot(Xh,Yh,'k--','LineWidth',1.);
	end
	end
end

% Calculate the limb, i.e. the edge of the planet disc on the sky
the = linspace(0,pi,200);

discrim = abs(cos(the).*cos(theobs))-(1-ecc^2)*abs(sin(the)*sin(theobs));
f=figure; plot(the,discrim), title('discrim limb'), pause, close(f)
the = the(discrim < 0);

phi_rel = acos(-cos(the).*cos(theobs)./((1-ecc^2)*sin(the)*sin(theobs)));
phi = phi_rel + phiobs;
f=figure; plot(the,phi_rel), pause, close(f)

phi1 =  phi_rel + phiobs;
phi2 = -phi_rel + phiobs;

if 1,
phi = [phi1,fliplr(phi2)];
the = [the,fliplr(the)];
else
phi = [phi1];
the = [the];
end
f=figure; plot(the,phi), title('limb the(phi)'), pause, close(f)

[r,x,y,z,xsky,ysky,zsky] = spherical2Sky(semiMaj_km,ecc, ...
                                    the,phi,theobs,phiobs,km2asec);
plotSky(theobs,phiobs,x,y,z,xsky,ysky,zsky)

if 0,
xsky(zsky<0) = NaN;
ysky(zsky<0) = NaN;
end
[min(zsky),max(zsky)]
% rotation
[xsky,ysky] = Rotate(alpha,xsky,ysky);
%plot(xsky, ysky, 'r-', -xsky, ysky, 'r-','LineWidth',1);
plot(xsky, ysky, 'r-', 'LineWidth',2);
pause

%return

% Repeat the above calculation for the locus of the day/night terminator
the = linspace(0,pi,200);

discrim = abs(cos(the).*cos(thesun))-(1-ecc^2)*abs(sin(the)*sin(thesun));
f=figure; plot(the,discrim), title('discrim terminator'), pause, close(f)
the = the(discrim < 0);

phi_rel = acos(-cos(the).*cos(thesun)./((1-ecc^2)*sin(the)*sin(thesun)));

% Retrieve the longitudes of the points on the 
% terminator, taking into account the difference between
% the longitude of the Earth (obs) and Sun

fprintf(1,'Sun is at relative longitude %.4f to sub-Earth point\n', ssedlon);

phi1 =  phi_rel + phiobs + ssedlon*pi/180.;
phi2 = -phi_rel + phiobs + ssedlon*pi/180.;

if 1,
phi = [phi1,fliplr(phi2)];
the = [the,fliplr(the)];
else
phi = [phi1];
the = [the];
end
f=figure; plot(the,phi), title('terminator the(phi)'), pause, close(f)

[r,x,y,z,xsky,ysky,zsky] = spherical2Sky(semiMaj_km,ecc, ...
                                         the,phi,theobs,phiobs,km2asec);

plotSky(theobs,phiobs,x,y,z,xsky,ysky,zsky)

if 0,
xsky(zsky*sign(ssedlon)<0) = NaN;
ysky(zsky*sign(ssedlon)<0) = NaN;
end
[min(zsky),max(zsky)]
% rotation
[xsky,ysky] = Rotate(alpha,xsky,ysky);
Xv = xsky; Xv(zsky*sign(grl2srl(ssedlon))<0) = NaN;
Yv = ysky; Yv(zsky*sign(grl2srl(ssedlon))<0) = NaN;
Xh = xsky; Xh(zsky*sign(grl2srl(ssedlon))>=0) = NaN;
Yh = ysky; Yh(zsky*sign(grl2srl(ssedlon))>=0) = NaN;
plot(Xv, Yv, 'g-','LineWidth',2);
plot(Xh, Yh, 'g--','LineWidth',2);


% limb and terminator
if 1,
f=figure;
r = semiMaj_km*km2asec;
[xll,yll,xld,yld,xtl,ytl,xtd,ytd,xc,yc] = ltc(r,ecc,se,ss);
fprintf(1,'cusp points x=%f,%f, y=%f,%f\n',xc',yc');
close(f);
f=figure;
[ll,ld,tl,td,cusp]=getLTC(r,ecc,se,ss);
fprintf(1,'cusp points x=%f,%f, y=%f,%f\n',cusp{1},cusp{2});
close(f);

% rotation
[xc,yc] = Rotate(alpha,cusp{1},cusp{2});
plot(xc,yc,'bx-','LineWidth',1,'MarkerSize',5);
end

hold off

function [r,x,y,z,xsky,ysky,varargout]=spherical2Sky(a,e,the,phi,theobs,phiobs,km2asec)

% ellipsoid with flattening along z 
r = a*sqrt(1-e^2)./sqrt(1-(e*sin(the)).^2);
if length(r)==1, % trick to make z same size as x and y
  r = r*ones(size(phi));
end

x = r .* sin(the) .* cos(phi);
y = r .* sin(the) .* sin(phi);
z = r .* cos(the);

Pm = [ sin(phiobs)            , cos(phiobs)            ,0          ;...
      -cos(phiobs)*cos(theobs),-sin(phiobs)*cos(theobs),sin(theobs);...
			 cos(phiobs)*sin(theobs), sin(phiobs)*sin(theobs),cos(theobs)];

A = viewmtx(90-phiobs*180/pi,90-theobs*180/pi);

% Project planetary position onto the plane of the sky
xsky = -x*sin(phiobs) + y*cos(phiobs);
ysky = (-x*cos(phiobs)-y*sin(phiobs))*cos(theobs)+z*sin(theobs);

% Transform to arc seconds 
xsky = km2asec*xsky;
ysky = km2asec*ysky;

if nargout > 6,
zsky = (x*cos(phiobs)+y*sin(phiobs))*sin(theobs)+z*cos(theobs);
varargout(1) = {zsky};
end

return

function [Vx,Vy] = Rotate(alpha,Ux,Uy)

% alpha is in deg
cosa = cosd(alpha);
sina = sind(alpha);

Vx = cosa*Ux - sina*Uy;
Vy = sina*Ux + cosa*Uy;

function plotSky(theobs,phiobs,x,y,z,xsky,ysky,zsky)

f=figure;

subplot(211),
plot3(x,y,z),
axis equal,
xlabel('x'),ylabel('y'),zlabel('z')
title(sprintf('theobs %.1f phiobs %.1f\n',180/pi*[theobs,phiobs]))

subplot(212),
ip = zsky>0;
plot3(xsky(ip),ysky(ip),zsky(ip),'-',...
      xsky(~ip),ysky(~ip),zsky(~ip),'--'),

%axis equal,
xlabel('x'),ylabel('y'),zlabel('z')
pause
close(f)

function lon=grl2srl(lon)
lon(lon>180) = lon(lon>180)-360;

