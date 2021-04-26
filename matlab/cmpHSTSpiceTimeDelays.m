function cmpHSTSpiceTimeDelays(filename)
% function cmpHSTSpiceTimeDelays(filename)

%
% $Id: cmpHSTSpiceTimeDelays.m,v 1.6 2021/04/26 13:27:59 patrick Exp $
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

params = getDefaultVOISEParams;
params.iFile = filename;
params.verbose = false;

% get HST image and parameters 
params = loadImage(params);

HST = params.HST;

planet = setUpSpice4Planet(HST);

% conversion factors
degPerRad = cspice_dpr;
radPerDeg = cspice_rpd;
radPerArcsec = radPerDeg/3600;
arcsecPerRad = degPerRad*3600;

% start and end exposure times
if ~isempty(HST.TDATEOBS) && ~isempty(HST.TTIMEOBS),
  % Convert string date to double precision 
  et = cspice_str2et([HST.TDATEOBS ' ' HST.TTIMEOBS]);
  if ~isempty(HST.EXPTIME),
    times    = et + [0, HST.EXPTIME]; 
  elseif isempty(HST.TEXPTIME),
    times    = et + [0, HST.TEXPTIME]; 
  end 
elseif ~isempty(HST.START_EPOCH) && ~isempty(HST.END_EPOCH)
  times = [cspice_str2et(HST.START_EPOCH), cspice_str2et(HST.END_EPOCH)];
end
format = 'C'; % Calendar format, UTC
%format = 'D'; % Day-of-Year format, UTC
%format = 'J'; % Julian Date format, UTC
%format = 'ISOC';  % ISO Calendar format, UTC
%format = 'ISOD';  % ISO Day-of-Year format, UTC
prec = 2; % number of decimal precision for fractional seconds computed
utcstr = cspice_et2utc(times,format,prec);
fprintf('Times             =   %s :   %s\n',utcstr(1,:),utcstr(2,:));

% get planet ra,dec in J2000 frame for difference aberration corrections
target   = planet.name;
frame    = 'J2000'; % Earth mean equator, dynamical equinox of J2000
observer = 'EARTH';

% aberration corrections for "reception" case (photons depart from the
% target's location at the light-time corrected epoch et-lt and *arrive* at the
% observer's location at 'et'
abcorr = { ...
  'NONE', ... % No correction
  'LT', ... % Correct for one-way light time
	'LT+S', ... % Correct for one-way light time and stellar aberration
	'CN', ... % Converged Newtonian light time correction
	'CN+S' ... % Converged Newtonian light time and stellar aberration
	};

range = zeros(length(times), length(abcorr));
lt = zeros(length(times), length(abcorr));
ra = zeros(length(times), length(abcorr));
dec = zeros(length(times), length(abcorr));

for i = 1:length(abcorr),
  [state, lt(:,i)] = cspice_spkezr(target, times, frame, abcorr{i}, observer);
  [range(:,i),ra(:,i),dec(:,i)] = cspice_recrad(state(1:3,:));
end

% rad to deg
ra = ra*degPerRad;
dec = dec*degPerRad;

%fprintf(1,'planet    ra, dec = %.4f,%.4f : %.4f,%.4f\n',[ra(:)';dec(:)']);
for i = 1:length(abcorr),
  fprintf(1,'pc ra, dec (deg)  = %12.6f,%12.6f : %12.6f,%12.6f (%s)\n',...
	[ra(:,i),dec(:,i)]', abcorr{i});
end

img = params.W;

% Mesh in pixels [1,nx]x[1,ny]
[Xj,Yj]=meshgrid(params.x+1,params.y+1);

% ref pixel
rpx = HST.CRPIX1;
rpy = HST.CRPIX2;
% ref ra/dec
rpra = HST.CRVAL1;
rpdec = HST.CRVAL2;

% pixel coordinates to world coordinates (ra/dec) (indices i)
[Xi,Yi] = getHSTpixel2radec(HST,Xj,Yj);

% planet world coordinates (ra/dec) to pixel coordinates
[pxc,pyc] = getHSTradec2pixel(HST,ra,dec);

for i = 1:length(abcorr),
  fprintf(1,'pc  x, y (pixel)  = %12.6f,%12.6f : %12.6f,%12.6f (%s)\n',...
	[pxc(:,i),pyc(:,i)]', abcorr{i});
end


close all

opts1 = {'fontsize',9,'fontweight','light','color','black'};
opts2 = {'fontsize',9,'fontweight','light','color','magenta'};

figure
if 1,
  [Xj,Yj] = getHSTabs2relPixels(HST,Xj,Yj);
	[rpx,rpy] = getHSTabs2relPixels(HST,rpx,rpy);
	[pxc,pyc] = getHSTabs2relPixels(HST,pxc,pyc);
	xlbl = 'pixels/ref. pixel [arcsec]';
	ylbl = 'pixels/ref. pixel [arcsec]';
else
	xlbl = 'pixels';
	ylbl = 'pixels';
end

pcolor(Xj,Yj,log10(abs(img))); shading flat;
hold on
plot(rpx,rpy,'ko','markersize',5);
plot(pxc(:,1:3),pyc(:,1:3),'-kx','markersize',5);
plot(pxc(:,4:5),pyc(:,4:5),':mo','markersize',5);
for i = 1:3,
  text(pxc(1,i)*1.05, pyc(1,i)*1.05,['S(' abcorr{i} ')'],opts1{:})
end
for i = 4:5,
  text(pxc(2,i)*1.05, pyc(2,i)*1.05,['E(' abcorr{i} ')'],opts2{:})
end
xlabel(xlbl)
ylabel(ylbl)
hold off
axis auto
axis equal

figure
if 1,
  % convert to arcsec and set relative ra/dec to ref pix
	[Xi,Yi] = getHSTabs2relRadec(HST,Xi,Yi);
	[rpra,rpdec] = getHSTabs2relRadec(HST,rpra,rpdec);
	[ra,dec] = getHSTabs2relRadec(HST,ra,dec);
	xlbl = 'ra/ref. pixel [arcsec]';
	ylbl = 'dec/ref. pixel [arcsec]';
else
  % absolute ra/dec in deg
	xlbl = 'ra [deg]';
	ylbl = 'dec [deg]';
end

pcolor(Xi,Yi,log10(abs(img))); shading flat;
hold on
axis auto
plot(rpra,rpdec,'ko','markersize',5)
plot(ra(:,1:3),dec(:,1:3),'-kx','markersize',5)
plot(ra(:,4:5),dec(:,4:5),':mo','markersize',5)
for i = 1:3,
  text(ra(1,i)*1.05, dec(1,i)*1.05,['S(' abcorr{i} ')'],opts1{:})
end
for i = 4:5,
  text(ra(2,i)*1.05, dec(2,i)*1.05,['E(' abcorr{i} ')'],opts2{:})
end
hold off
axis auto
axis equal
xlabel(xlbl);
ylabel(ylbl);


%  particularly in MATLAB due to data persistence.
cspice_kclear

