function plotPlanetGrid(planet,params,pc,epoch,CML,psi,orientat,PIXSIZE)
% function plotPlanetGrid(planet,params,pc,epoch,CML,psi,orientat,PIXSIZE)

%
% $Id: plotPlanetGrid.m,v 1.2 2015/09/25 17:55:09 patrick Exp $
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


function drawPlanetGrid(planet,epoch,ss,se,orientat,dlat,dlon,semiMaj_km,ecc)

% Get difference between (right-handed) sub-Earth and sub-solar longitudes.
% Although these two quantities will vary on time scale of the planet's
% rotation, their *difference* varies much more slowly. 
% Note that a *positive* value means that the Sun central meridian is
% situated ahead of the Earth's (in the direction of planetary rotation)
ssedlon = ss.lon - se.lon;

% Sub-Earth colatitude in radians
theobs = pi/2 - se.lat*pi/180;
% Central meridian longitude in a right-handed system
phiobs = 2*pi - se.CML*pi/180;

% Sub-solar colatitude in radians
thesun = pi/2 - ss.lat*pi/180;

% Scaling factor to convert from km to arcsec on the plane of 
% the Earth observer's sky
km2asec = (1/se.distkm)*(180/pi)*3600;

if 0
mfig = gcf;
figure
r = semiMaj_km*km2asec;
[xll,yll,xld,yld,xtl,ytl,xtd,ytd,xc,yc] = ltc(r,ecc,se,ss);
fprintf(1,'cusp points x=%f,%f, y=%f,%f\n',xc',yc');
[ll,ld,tl,td,cusp]=getLTC(r,ecc,se,ss);
fprintf(1,'cusp points x=%f,%f, y=%f,%f\n',cusp{1},cusp{2});
figure(mfig)
end

% Calculate the viewing angle
lineOfSight = [sin(theobs)*cos(phiobs);sin(theobs)*sin(phiobs);cos(theobs)];

% http://www.imcce.fr/en/ephemerides/formulaire/form_ephephys.php
% gives 254.31 deg 
alpha = se.psi-orientat;
fprintf(1,'psi %f orientat %f alpha %f\n',se.psi,orientat,alpha);

hold on

% Grid curves of constant latitude
for the = ([dlat:dlat:180-dlat])*pi/180,
%for the = ([90+se.lat])*pi/180,

	phi = linspace(0,2*pi,fix(150*sin(the))); 

	[r,x,y,z,xsky,ysky,zsky] = spherical2Sky(semiMaj_km,ecc, ...
	                                    the,phi,theobs,phiobs,km2asec);

  if 0,
  %cos_vang = lineOfSight(1)*x./r+lineOfSight(2)*y./r+lineOfSight(3)*z./r;
  cos_vang = lineOfSight(1)*x+lineOfSight(2)*y+lineOfSight(3)*z;
	flag_vang = (cos_vang > 0);
	else
	flag_vang = (zsky > 0);
	end
	if 0
	Xv = xsky(flag_vang);
	Yv = ysky(flag_vang);
	[Xv,I] = sort(Xv);
	Yv = Yv(I);
	else
	Xv = xsky; Xv(~flag_vang) = NaN;
	Yv = ysky; Yv(~flag_vang) = NaN;
	end
	% rotation
	[Xv,Yv] = Rotate(alpha,Xv,Yv);
  plot(Xv,Yv,'k-','LineWidth',.5);
	if 0, % not visible line
	if 0,
	Xh = xsky(~flag_vang);
	Yh = ysky(~flag_vang);
	[Xh,I] = sort(Xh);
	Yh = Yh(I);
	else
	Xh = xsky; Xh(flag_vang) = NaN;
	Yh = ysky; Yh(flag_vang) = NaN;
	end
  [Xh,Yh] = Rotate(alpha,Xh,Yh);
  plot(Xh,Yh,'k--','LineWidth',.5);
	end
end

% Grid curves of constant longitude
for phi = [0:dlon:360-dlon]*pi/180,

  the = linspace(0,pi,50);

	[r,x,y,z,xsky,ysky,zsky] = spherical2Sky(semiMaj_km,ecc, ...
	                                    the,phi,theobs,phiobs,km2asec);

  if 0,
  cos_vang = lineOfSight(1)*x./r+lineOfSight(2)*y./r+lineOfSight(3)*z./r;
  cos_vang = lineOfSight(1)*x+lineOfSight(2)*y+lineOfSight(3)*z;
  flag_vang = (cos_vang > 0);
  else
	flag_vang = (zsky > 0);
	end
	if 0, 
	Xv = xsky(flag_vang);
	Yv = ysky(flag_vang);
	[Yv,I] = sort(Yv);
	Xv = Xv(I);
	else
	Xv = xsky; Xv(~flag_vang) = NaN;
	Yv = ysky; Yv(~flag_vang) = NaN;
	end
	% rotation
	[Xv,Yv] = Rotate(alpha,Xv,Yv);
	if phi>0,
  plot(Xv,Yv,'k-','LineWidth',.5);
	else, % phi==0,
	plot(Xv,Yv,'c-','LineWidth',1.0);
	end
	if 0,
	if 0,
	Xh = xsky(~flag_vang);
	Yh = ysky(~flag_vang);
	[Yh,I] = sort(Yh);
	Xh = Xh(I);
	else
	Xh = xsky; Xh(flag_vang) = NaN;
	Yh = ysky; Yh(flag_vang) = NaN;
	end
  [Xh,Yh] = Rotate(alpha,Xh,Yh);
  plot(Xh,Yh,'k--','LineWidth',.5);
	end
end

% Calculate the limb, i.e. the edge of the planet disc on the sky
the = linspace(0,pi,200);

discrim = abs(cos(the).*cos(theobs))-(1-ecc^2)*abs(sin(the)*sin(theobs));
f=figure; plot(the,discrim), pause, close(f)
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
f=figure; plot(the,phi), pause, close(f)

[r,x,y,z,xsky,ysky,zsky] = spherical2Sky(semiMaj_km,ecc, ...
                                    the,phi,theobs,phiobs,km2asec);

if 1,
xsky(zsky<0) = NaN;
ysky(zsky<0) = NaN;
end

f=figure;
subplot(211), plot3(x,y,z),
subplot(212), plot3(xsky,ysky,zsky),
pause, close(f)

% rotation
[xsky,ysky] = Rotate(alpha,xsky,ysky);
%plot(xsky, ysky, 'r-', -xsky, ysky, 'r-','LineWidth',1);
plot(xsky, ysky, 'r-', 'LineWidth',1);
pause

%return

% Repeat the above calculation for the locus of the day/night terminator
the = linspace(0,pi,200);

discrim = abs(cos(the).*cos(thesun))-(1-ecc^2)*abs(sin(the)*sin(thesun));
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
f=figure; plot(the,phi), pause, close(f)

[r,x,y,z,xsky,ysky,zsky] = spherical2Sky(semiMaj_km,ecc, ...
                                         the,phi,theobs,phiobs,km2asec);

if 0,
flag_zsky = (zsky*sign(ssedlon) > 0);
xsky = xsky(flag_zsky);
ysky = ysky(flag_zsky);
end

if 0,
% trick to avoid "jumps" in theta
dt=diff(atan2(ysky,xsky));
ii=find(abs(dt)>2*min(abs(dt)));
xsky(ii) = NaN;
ysky(ii) = NaN;
end

if 1,
xsky(zsky*sign(ssedlon)<0) = NaN;
ysky(zsky*sign(ssedlon)<0) = NaN;
end

f=figure;
subplot(211), plot3(x,y,z),
subplot(212), plot3(xsky,ysky,zsky),
pause, close(f)

% rotation
[xsky,ysky] = Rotate(alpha,xsky,ysky);
plot(xsky, ysky, 'g-','LineWidth',1);

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


% Project planetary position onto the plane of the sky
xsky = -x*sin(phiobs) + y*cos(phiobs);
ysky = (-x*cos(phiobs)-y*sin(phiobs))*cos(theobs)+z*sin(theobs);

% Transform to arc seconds 
xsky = km2asec*xsky;
ysky = km2asec*ysky;

if 0
f=figure;
fprintf(1,'theobs %f phiobs %f\n',180/pi*[theobs,phiobs])
subplot(211),plot3(x,y,z),axis equal,xlabel('x'),ylabel('y'),zlabel('z')
zsky = (x*cos(phiobs)+y*sin(phiobs))*sin(theobs)+z*cos(theobs);
zsky = km2asec*zsky;
ip = zsky>0;
subplot(212),plot3(xsky(ip),ysky(ip),zsky(ip),'-',...
                   xsky(~ip),ysky(~ip),zsky(~ip),'--'),
axis equal, xlabel('x'),ylabel('y'),zlabel('z')
pause
close(f)
end


if nargout > 6,
zsky = (x*cos(phiobs)+y*sin(phiobs))*sin(theobs)+z*cos(theobs);
varargout(1) = {zsky};
end

return

function [Vx,Vy] = Rotate(alpha,Ux,Uy)

cosa = cosd(alpha);
sina = sind(alpha);

Vx = cosa*Ux - sina*Uy;
Vy = sina*Ux + cosa*Uy;


