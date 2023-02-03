function testPolProjAPISproc(filename)
% function testPolProjAPISproc(filename)
%
% polproj('../share/input/j9rlb0fxq_proc.fits')

%
% $Id: testPolProjAPISproc.m,v 1.1 2023/02/03 16:48:28 patrick Exp $
%
% Copyright (c) 2021 Patrick Guio <patrick.guio@gmail.com>
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

close all

verbose = true;

params = getDefaultVOISEParams;
params.iFile = filename;
params.verbose = verbose;

params = loadImage(params);

%close all

binlat = 1;
dawnmax = 13;
duskmin = 11;
np = 3;
params = minnaert(params,binlat,dawnmax,duskmin,np);
params = li(params,binlat,dawnmax,duskmin);

%subplot(211)
%imagesc(params.x,params.y,params.W)
%axis xy

figure

%subplot(212)
%pcolor(params.x,params.y,params.W)
%shading flat

if 1,
latgrid = params.ext.lat1b;
ltgrid = params.ext.lt1b;
else
latgrid = params.ext.lat300km;
ltgrid = params.ext.lt300km;
end

maskPlanet = find(latgrid~=-100);
maskNotPlanet = find(latgrid==-100);

maskPlanet = find(latgrid<=-45 & latgrid~=-100);
maskNotPlanet = find(latgrid>-45 | latgrid==-100);

W1 = params.W-params.minnaert.airglow;
W2 = params.W-params.li.airglow;
clim = [0,max([W1(:);W2(:)])];
fprintf(1,'clim = %f,%f\n',clim);

W1(maskNotPlanet) = NaN;
W2(maskNotPlanet) = NaN;

xlim = [60,1020];
ylim=[560,760];

subplot(311),
pcolor(params.x,params.y,W1)
shading flat
set(gca,'xlim',xlim,'ylim',ylim);
set(gca,'clim',clim);
colorbar
title(['Minnaert np=' num2str(np)])

subplot(312),
pcolor(params.x,params.y,W2)
shading flat
set(gca,'xlim',xlim,'ylim',ylim);
set(gca,'clim',clim);
colorbar
title('Li')

subplot(313),
pcolor(params.x,params.y,W1-W2)
shading flat
set(gca,'xlim',xlim,'ylim',ylim);
[~,dmin,dmax,~] = isoutlier(W1(:)-W2(:),'median');
dlim = 3*[dmin,dmax];
set(gca,'clim',dlim);
colorbar
title('Minnaert-Li')

orient tall
print('-dpdf',replace(filename,'.fits','_1.pdf'));

pause

fprintf(1,'latitude extent  = %+8.1f, %+8.1f deg\n', ...
        [min(latgrid(maskPlanet)),max(latgrid(maskPlanet))])
% local time = 12 longitude =   0
% local time = 06 longitude =  90 dawn
% local time = 18 longitude = -90/270 dusk
% local time = 00 longitude = 180
longrid = lt2lon(ltgrid);
fprintf(1,'longitude extent = %+8.1f, %+8.1f deg\n', ...
        [min(longrid(maskPlanet)),max(longrid(maskPlanet))])

[x,y,z] = latlon2xyz(latgrid(maskPlanet),longrid(maskPlanet));
c1 = W1(maskPlanet);
c2 = W2(maskPlanet);

% griddedInterpolant({x1g,x2g,...,xng}
opts = {'linear','none'};
F1 = scatteredInterpolant([latgrid(maskPlanet),longrid(maskPlanet)],c1,opts{:});

F2 = scatteredInterpolant([latgrid(maskPlanet),longrid(maskPlanet)],c2,opts{:});

%[xgrid,ygrid,zgrid]=latlon2xyz(latgrid(maskPlanet),longrid(maskPlanet));
%opts = {'linear','none'};
%Fxyz = scatteredInterpolant(xgrid,ygrid,zgrid,c,opts{:});

xi = linspace(-0.8,0.8,200);
yi = linspace(-0.8,0.8,200);
e = 0.3543;
[X,Y] = meshgrid(xi,yi);
Z = -sqrt(1-e^2)*sqrt(1-X.^2-Y.^2);
ii = find(X.^2-Y.^2>=1);
ii = find(imag(Z)~=0);
X(ii) = 1; Y(ii) = 0; Z(ii) = 0;
[LAT,LON] = xyz2latlon(X,Y,Z);
C1 = F1(LAT,LON);
C2 = F2(LAT,LON);

subplot(311),
pcolor(X,Y,C1),
shading flat
set(gca,'xlim',[-.8,.8],'ylim',[-.8,.8])
axis square
colorbar
set(gca,'clim',clim);
%pause

subplot(312),
pcolor(X,Y,C2),
shading flat
set(gca,'xlim',[-.8,.8],'ylim',[-.8,.8])
axis square
colorbar
set(gca,'clim',clim);

subplot(313),
pcolor(X,Y,C1-C2),
shading flat
set(gca,'xlim',[-.8,.8],'ylim',[-.8,.8])
axis square
colorbar
set(gca,'clim',dlim);

orient tall
print('-dpdf',replace(filename,'.fits','_2.pdf'));


figure

subplot(311)
polarproj(x,y,c1,clim,['Minnaert np=' num2str(np)])

subplot(312)
polarproj(x,y,c2,clim,'Li')

subplot(313)
polarproj(x,y,c1-c2,dlim,'Minneart-Li')

orient tall
print('-dpdf',replace(filename,'.fits','_3.pdf'));

if 0,
[x,y,z] = latlon2xyz([-30,-45],60)
[lat,lon] = xyz2latlon(x,y,z)
end

return

x = linspace(-0.2,0.0,50);
y = linspace(0.0,0.4,50);
[X,Y] = meshgrid(x,y);
[LAT,LON] = xyz2latlon(X,Y,0);
imagesc(x,y,LAT),colorbar,pause
imagesc(x,y,LON),colorbar,pause
for i=1:50,
   for j=1:50,
     dist = sqrt((latgrid-LAT(i,j)).^2+(longrid-LON(i,j)).^2);
		 [mn,id] = min(dist(:))
		 imagesc(x,y,dist)
		 colorbar
		 drawnow
		 %pause
		end
end

function y = model(b,x)
% Jupiter
e = 0.3543;
re = 71492; % km
rp = 66854; % km
lat = x(:,:,1);
lon = x(:,:,2);
x = b(1)*sind(lon).*cosd(lat);
y = b(1)*cosd(lon).*cosd(lat);
z = b(1)*sqrt(1-e^2)*sind(lat);

function lon = lt2lon(lt)
lon = 180-lt/24*360;

function [x,y,z] = latlon2xyz(lat,lon)
% Jupiter
e = 0.3543;
re = 71492; % km
rp = 66854; % km
x = sind(lon).*cosd(lat);
y = cosd(lon).*cosd(lat);
z = sqrt(1-e^2)*sind(lat);

function [lat,lon] = xyz2latlon(x,y,z)
% Jupiter
e = 0.3543;
re = 71492; % km
rp = 66854; % km
%lat = -acosd(sqrt(x.^2+y.^2));
lat = atan2d(z/sqrt(1-e^2),sqrt(x.^2+y.^2));
lon = atan2d(x,y);

function polarproj(x,y,c,clim,model)

scatter(x,y,[],c,'filled')
hold on
lat = -[20:10:80,85];
lat = -[45:10:85];
lon = linspace(0,360,100);
for i=1:length(lat),
  [x,y,z] = latlon2xyz(lat(i),lon);
  plot(x,y,'k:')
end
lat = -linspace(20,85,100);
lat = -linspace(45,85,100);
lon = [0:15:360];
for i=1:length(lon),
  [x,y,z] = latlon2xyz(lat,lon(i));
  plot(x,y,'k:')
end
[x,y,z] = latlon2xyz([-90,-45],0);
noon = plot(x,y,'r-');
[x,y,z] = latlon2xyz([-90,-45],90);
dawn = plot(x,y,'c-');
[x,y,z] = latlon2xyz([-90,-45],-90);
dusk = plot(x,y,'m-');
[x,y,z] = latlon2xyz([-90,-45],180);
midnight = plot(x,y,'k-');
title(model)
legend([noon,dawn,dusk,midnight],{'noon','dawn','dusk','midnight'},...
       'location','eastoutside')
axis square
hold off
colorbar
set(gca,'clim',clim);

