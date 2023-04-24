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

close all

binlat = 1;
dawnmax = 13;
duskmin = 11;
np = 3;
params = minnaert(params,binlat,dawnmax,duskmin,np);
params = li(params,binlat,dawnmax,duskmin);


if 0
figure
subplot(211)
imagesc(params.x,params.y,params.W)
axis xy

subplot(212)
pcolor(params.x,params.y,params.W)
shading flat
end

if 1,
latgrid = params.ext.lat1b;
ltgrid = params.ext.lt1b;
else
latgrid = params.ext.lat300km;
ltgrid = params.ext.lt300km;
end



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
% Correction value added to cos(diff) * cos(spec) + corr
% (Different values may work better)
corr = 0.5; 
method = 'phong'; % Can be phong, lambertian or none (specularity model)
%method = 'lambertian';
% Creates reference model for specularity and diffusion when looking at the
% planet from the sub-earth point direction.
adj = intensityAdjustment(x, y, z, params.HST.APIS.SUBELAT, ...
    params.HST.APIS.SUBELON, params.HST.APIS.SUBSLAT, ...
    params.HST.APIS.SUBSLON, corr, method);
reorient = intensityAdjustment(x, y, z, -90, params.HST.APIS.SUBELON,...
    params.HST.APIS.SUBSLAT,params.HST.APIS.SUBSLON, 0, 'none');
% First two lines are dividing original images by the reference model such
% that the images' intensities are adjusted
c1 = c1 ./ adj;
c2 = c2 ./ adj;
% Multiplies with an assumed lambertian scatter for the observer at the
% south pole (Highest values directly at the pole.
c1 = c1 .* reorient; 
c2 = c2 .* reorient;
end

opts = {'linear', 'none'}; % Nearest instead of linear can produce slightly different result

% Few points around the pole, interpolation somewhat fixes this.
F1 = scatteredInterpolant(x, y, c1, opts{:});
F2 = scatteredInterpolant(x, y, c2, opts{:});

xi = linspace(min(x), max(x), 400); % Seems to be trivial difference above 400
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
dlim = conversionFactor(dlim);

pause

showLogarithmic = 0; % Set to 1 to display in logarithmic scale
figure
subplot(311)
polarproj(X, Y, C1, clim, 'Minnaert', showLogarithmic) 
subplot(312)
polarproj(X, Y, C2, clim, 'Li', showLogarithmic)
subplot(313)
polarproj(X, Y, C1 - C2, dlim, 'Minnaert - Li', showLogarithmic)

orient tall
print('-dpdf',replace(filename,'.fits','_2.pdf'));

pause

% Alternate plotting using polarscatter instead of scatter
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
re = 71492; % km
rp = 66854; % km
e = sqrt(1 - rp^2 / re^2);
lat = x(:,:,1);
lon = x(:,:,2);
x = b(1)*sind(lon).*cosd(lat);
y = b(1)*cosd(lon).*cosd(lat);
z = b(1)*sqrt(1-e^2)*sind(lat);

function lon = lt2lon(lt)
lon = 180-lt/24*360;

function [x,y,z] = latlon2xyz(lat,lon)
% Jupiter
re = 71492; % km
rp = 66854; % km
e = sqrt(1 - rp^2 / re^2);
% Eq. 4 in planetproj
x = re * cosd(lon) .* cosd(lat);
y = re * sind(lon) .* cosd(lat);
z = re * sqrt(1 - e^2) * sind(lat);

function [lat,lon] = xyz2latlon(x,y,z)
% Jupiter
re = 71492; % km
rp = 66854; % km
e = sqrt(1 - rp^2 / re^2);
%lat = -acosd(sqrt(x.^2+y.^2));
lat = atan2d(z/sqrt(1-e^2),sqrt(x.^2+y.^2));
lon = atan2d(x,y);

function latCorrected = planetodetic(lat)
% Takes into account local gravity
% Approximation using Clairaut's equation
% Eq. 6 and 7 in planetproj 
re = 71492; % km
rp = 66854; % km
f = (re - rp) / re;
latCorrected = lat + atand(f * sind(2 .* lat) ./ (1 - f .* sind(lat).^2));

function polarproj(x, y, c, clim, model, logScale)
% Plots values using cartesian coordinates
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
% Plotgrid made own function for better readability
% Plots gridlines
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
% Alternate display using polarscatter which uses polar coordinates
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

function adjustment = intensityAdjustment(x, y, z, subelat, subelon, ...
    subslat, subslon, corr, method)
% Adjusts the intensity taking into account the viewing angle and incidence
% angle. The specularity can be corrected for as well with two different
% models. The normal lambertian scatter model for the sun hitting the
% planet's surface and a phong model, which is dependent on a shininess
% exponent n.

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

% Scattered light from the surface to the observer
% Changed variable name so View vector in phong can be calculated correctly
subEarthPointNorm = subEarthPoint ./ norm(subEarthPoint);
dotProdEm = Nhat * subEarthPointNorm'; % Emission
emAngle = acosd(dotProdEm);

switch method % Method for specularity
    case 'phong'
        % https://dl.acm.org/doi/abs/10.1145/360825.360839
        % https://dicklyon.com/tech/Graphics/Phong_TR-Lyon.pdf
   
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
        
        n = 2; % Shininess exponent (Different values may work better)

        % Diffuse reflection coefficient
        Kdiff = dot(Lhat, Nhat, 2); % Seems to get an unexpected result
        
        % Specular reflection coefficient
        Kspec = dot(Rhat, Vhat, 2) .^ n;
        
        % Total reflection
        Ks = 1 ; % Specular reflection constant (assumed 1)
        Kd = 1 ; % Specular diffusion constant (assumed 1)
        
        Ktot = Kd .* Kdiff + Ks .* Kspec;
                
        % Total Illumination
        % Corr is value added to avoid large values near terminator and limb.
        % 1 avoids values larger than original, lower values makes the brightness
        % proportionally higher than the original maximum value, which lowers the
        % intensity in the rest of the image relatively if clim is not applied to
        % the plots. A too low value will make too much of the original image's
        % intensity be above the maximum limit.

        adjustment = cosd(emAngle) .* Ktot + corr;

    case 'lambertian'        
        subSolarPoint = subSolarPoint ./ norm(subSolarPoint);
        dotProdInc = Nhat * subSolarPoint'; % Incidence
        incAngle = acosd(dotProdInc);
        
        % Adds factor to avoid extremely high values near
        % terminator and limb.
        adjustment = cosd(emAngle) .* cosd(incAngle) + corr; 

    case 'none'
        adjustment = cosd(emAngle) + corr; % No specularity

    otherwise
        disp('method must be either phong or lambertian, or none for no specularity')
end

function c = conversionFactor(c)
c = c / (1.473 * 10^-3); % Convert to kR