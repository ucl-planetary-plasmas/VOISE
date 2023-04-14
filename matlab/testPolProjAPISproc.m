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

subplot(211)
imagesc(params.x, params.y, params.W)
axis xy

subplot(212)
pcolor(params.x, params.y, params.W)
shading flat

pause
close all

binlat = 1;
dawnmax = 13;
duskmin = 11;
np = 3;
params = minnaert(params, binlat, dawnmax, duskmin, np);
params = li(params, binlat, dawnmax, duskmin);

close all

if 1
latgrid = params.ext.lat1b;
ltgrid = params.ext.lt1b;
else
latgrid = params.ext.lat300km;
ltgrid = params.ext.lt300km;
end

maskPlanet = find(latgrid <= -45 & latgrid ~= -100);
maskNotPlanet = find(latgrid > -45 | latgrid == -100);

W1 = params.W - params.minnaert.airglow;
W2 = params.W - params.li.airglow;

W1(maskNotPlanet) = NaN;
W2(maskNotPlanet) = NaN;

clim = [0, max([W1(:); W2(:)])];
fprintf(1, 'clim = %f,%f\n', clim);

xlim = [60,1020];
ylim = [560,760];

[~, dmin, dmax, ~] = isoutlier(W1(:) - W2(:), 'median');
dlim = 3 * [dmin, dmax];
if 0
subplot(311)
pcolor(params.x, params.y, W1)
shading flat
set(gca, 'xlim', xlim, 'ylim', ylim);
set(gca, 'clim', clim);
colorbar
title(['Minnaert np=' num2str(np)])

subplot(312)
pcolor(params.x, params.y, W2)
shading flat
set(gca, 'xlim', xlim, 'ylim', ylim);
set(gca, 'clim', clim);
colorbar
title('Li')

subplot(313)
pcolor(params.x, params.y, W1 - W2)
shading flat
set(gca, 'xlim', xlim, 'ylim', ylim);
set(gca, 'clim', dlim);
colorbar
title('Minnaert-Li')

orient tall
print('-dpdf',replace(filename,'.fits','_1.pdf'));

pause
close all
end

fprintf(1,'latitude extent  = %+8.1f, %+8.1f deg\n', ...
        [min(latgrid(maskPlanet)), max(latgrid(maskPlanet))])
% local time = 12 longitude =   0
% local time = 06 longitude =  90 dawn
% local time = 18 longitude = -90/270 dusk
% local time = 00 longitude = 180
longrid = lt2lon(ltgrid);
fprintf(1,'longitude extent = %+8.1f, %+8.1f deg\n', ...
        [min(longrid(maskPlanet)), max(longrid(maskPlanet))])

c1 = W1(maskPlanet);
c2 = W2(maskPlanet);

% Corrected for local gravity field
latgrid = planetodetic(latgrid(maskPlanet));
[x, y, z] = latlon2xyz(latgrid, longrid(maskPlanet));

% Increases intensity as the angle between observer and surface normal
% increases. When plotting, clim limits max intensity. Set to 0 for
% original intensities in plot.
adjust_intensity = 1;

if adjust_intensity
%subELAT = params.HST.APIS.SUBELAT;
subELAT = -90; % Viewed from the south-pole
adj = intensityAdjustment(x, y, z, subELAT, params.HST.APIS.SUBELON,...
    params.HST.APIS.SUBSLAT,params.HST.APIS.SUBSLON);
c1 = c1 .* adj;
c2 = c2 .* adj;
end

opts = {'linear', 'none'};

%c2_corr = params.ext.oza1b(maskPlanet);

% Few points around the pole, interpolation somewhat fixes this.
F1 = scatteredInterpolant(x, y, c1, opts{:});
F2 = scatteredInterpolant(x, y, c2, opts{:});

xi = linspace(min(x), max(x), 400);
yi = linspace(min(y), max(y), 400);
[X, Y] = meshgrid(xi, yi);

% rho calculated for limits in polar plot
[~, rho] = cart2pol(x, y);
[THETA, RHO] = cart2pol(X, Y);

C1 = F1(X, Y);
C2 = F2(X, Y); 

%Convert from electrons S^-1 to kR
C1 = conversionFactor(C1);
C2 = conversionFactor(C2);
clim = conversionFactor(clim);
%dlim = conversionFactor(dlim);

showLogarithmic = 1;
figure
subplot(311)
polarproj(X, Y, C1, clim, 'Minnaert', showLogarithmic) 
subplot(312)
polarproj(X, Y, C2, clim, 'Li', showLogarithmic)
subplot(313)
polarproj(X, Y, C1 - C2, clim, 'Minnaert - Li', showLogarithmic)

orient tall
print('-dpdf',replace(filename,'.fits','_2.pdf'));

% Alternate plotting using polarscatter instead of scatter
% Not logarithmic for comparison
showLogarithmic = 0;
figure
subplot(311)
polarprojAlt(THETA, RHO, C1, clim, 'Minnaert', rho, showLogarithmic)
subplot(312)
polarprojAlt(THETA, RHO, C2, clim, 'Li', rho, showLogarithmic)
subplot(313)
polarprojAlt(THETA, RHO, C1 - C2, dlim, 'Minnaert - Li', rho, showLogarithmic)

orient tall
print('-dpdf',replace(filename,'.fits','_3.pdf'));

return

x = linspace(-0.2, 0.0, 50);
y = linspace(0.0, 0.4, 50);
[X, Y] = meshgrid(x, y);
[LAT, LON] = xyz2latlon(X, Y, 0);
imagesc(x,y,LAT), colorbar, pause
imagesc(x,y,LON), colorbar, pause
for i = 1: 50
   for j = 1: 50
     dist = sqrt((latgrid - LAT(i, j)).^2 + (longrid - LON(i, j)).^2);
		 [mn, id] = min(dist(:));
		 imagesc(x, y, dist)
		 colorbar
		 drawnow
		 %pause
		end
end

function c = conversionFactor(c)
c = c / (1.473 * 10^-3); % Convert to kR

function lon = lt2lon(lt)
lon = 180-lt/24*360;

function [x,y,z] = latlon2xyz(lat, lon)
% Jupiter
re = 71492; % km
rp = 66854; % km

e = sqrt(1 - rp^2 / re^2);
% Eq. 4 in planetproj
x = re * cosd(lon) .* cosd(lat);
y = re * sind(lon) .* cosd(lat);
z = re * sqrt(1 - e^2) * sind(lat);

function latCorrected = planetodetic(lat)
% Takes into account local gravity
% Approximation using Clairaut's equation
% Eq. 6 and 7 in planetproj 
re = 71492; % km
rp = 66854; % km
f = (re - rp) / re;
latCorrected = lat + atand(f * sind(2 .* lat) ./ (1 - f .* sind(lat).^2));

function [lat,lon] = xyz2latlon(x, y, z)
% Jupiter
re = 71492; % km
rp = 66854; % km
e = sqrt(1 - rp^2 / re^2);
lat = atan2d(z / sqrt(1 - e^2), sqrt(x.^2 + y.^2));
lon = atan2d(x, y);

function polarproj(x, y, c, clim, model, logScale)
scatter(x(:), y(:), [], c(:), 'filled')
hold on
plotgrid
hold off
axis square
colorbar
ylabel(colorbar, 'kR H_2')
set(gca,'clim',clim);
if logScale
set(gca,'ColorScale','log')
end
title(model)
view(90, -90)

function plotgrid
lat = -(45: 10: 85);
lat = planetodetic(lat);
lon = linspace(0, 360, 100);
for i = 1: length(lat)
  [x, y, ~] = latlon2xyz(lat(i), lon);
  plot(x, y, 'k:')
end
lat = -linspace(45, 85, 100);
lat = planetodetic(lat);
lon = 0: 15: 360;
for i = 1: length(lon)
  [x, y, ~] = latlon2xyz(lat, lon(i));
  plot(x, y, 'k:')
end
[x, y, ~] = latlon2xyz([-90, -45], 0);
noon = plot(x, y, 'r-');
[x, y, ~] = latlon2xyz([-90, -45], 90);
dawn = plot(x, y, 'c-');
[x, y, ~] = latlon2xyz([-90, -45], -90);
dusk = plot(x,y,'m-');
[x, y, ~] = latlon2xyz([-90, -45], 180);
midnight = plot(x, y, 'k-');
legend([noon, dawn, dusk, midnight], {'noon', 'dawn', 'dusk', 'midnight'},...
       'location', 'eastoutside')

function polarprojAlt(theta, rho, c, clim, model, rlimit, logScale)
polarscatter(theta(:), rho(:), [], c(:), 'filled')
hold on

rmax = max(rlimit); 
noon = polarplot([0, 0], [0, rmax], Color="r");
dawn = polarplot([pi / 2, pi / 2], [0, rmax], Color="c");
dusk = polarplot([3 / 2 * pi, 3 / 2 * pi], [0, rmax], Color="m");
midnight = polarplot([pi, pi], [0, rmax], Color="w");
legend([noon, dawn, dusk, midnight], {'noon', 'dawn', 'dusk', 'midnight'})

ax = gca;
ax.Layer = 'top';
ax.GridLineStyle = ':';
ax.Color = [0, 0, 0];
ax.GridColor = [1, 1, 1];
ax.GridAlpha = 0.5;
ax.ThetaZeroLocation = 'top';
ax.ThetaDir = 'clockwise';
ax.RTick = linspace(0, rmax, 6);
ax.CLim = clim;

rlim([0, rmax]);
rticklabels([]);
thetaticks(0:15:360);
thetaticklabels({0, '', '', '', '', '',...
    90, '', '', '', '', '',...
    180, '', '', '', '', '',...
    270, '', '', '', '', ''})

title(model)
colorbar
ylabel(colorbar, 'kR H_2')
if logScale
    set(gca,'ColorScale','log')
end

function adjustment = intensityAdjustment(x, y, z, subelat, subelon, subslat, subslon)
% https://en.wikipedia.org/wiki/Specular_highlight
% https://en.wikipedia.org/wiki/Phong_reflection_model
% https://dl.acm.org/doi/abs/10.1145/360825.360839
% https://dicklyon.com/tech/Graphics/Phong_TR-Lyon.pdf

n = 2; % Shininess exponent (Different values may work better)

% Jupiter
re = 71492; % km
rp = 66854; % km
e = sqrt(1 - rp^2 / re^2);

% Direction of normal surface
mag = 1 ./ sqrt(x.^2 + y.^2 + z.^2 / (1 - e^2)^2);
Nhat = mag .* [x, y, z ./ (1 - e^2)];

[earthx, earthy, earthz] = latlon2xyz(subelat, subelon);
subEarthPoint = [earthx, earthy, earthz];

[sunx, suny, sunz] = latlon2xyz(subslat, subslon);
subSolarPoint = [sunx, suny, sunz];

% View vector
% https://se.mathworks.com/help/matlab/ref/vecnorm.html
V = subEarthPoint - [x, y, z];
Vhat = V ./ vecnorm(V, 2, 2);

% Reflection vector
L = [x, y, z] - subSolarPoint;
Lhat = L ./ vecnorm(L, 2, 2);

% Direction Vector
% https://se.mathworks.com/help/matlab/ref/dot.html
R = 2 * dot(Nhat, Lhat, 2) .* Nhat - Lhat;
Rhat = R ./ vecnorm(R, 2, 2);

% Diffuse reflection coefficient
Kdiff = dot(Lhat, Nhat, 2); % Seems to get an unexpected result

% Specular reflection coefficient
Kspec = dot(Rhat, Vhat, 2) .^ n;

% Total reflection
Ks = 1 ; % Specular reflection constant (assumed 1)
Kd = 1 ; % Specular diffusion constant (assumed 1)

Ktotal = Kd .* Kdiff + Ks .* Kspec;

% Scattered light from the aurora
subEarthPoint = subEarthPoint ./ norm(subEarthPoint);
dotProd = Nhat * subEarthPoint';
theta = acosd(dotProd);

lambertianScatter = cosd(theta);

% Total Illumination
adjustment = lambertianScatter + Ktotal;
