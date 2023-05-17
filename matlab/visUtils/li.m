function params = li(params,binlat,dawnmax,duskmin)
% function params = li(params,binlat,dawnmax,duskmin)
%
% Compute a airglow model as described in
% Dayglow Removal from FUV Auroral Images, Li et al, 
% IGARSS, IEEE International Geoscience and Remote Sensing Symposium, 2004
%
%
% params : a VOISE parameter structure with image as created by
%          functions getDefaultVOISEParams and loadImage.
% binlat :  latitude bin size in degrees (default 1 deg)
% dawnmax: maximum local time for dawn sector fit (default 13h)
% duskmin: minimum local time for dusk sector fit (default 11h)

%
% $Id: li.m,v 1.2 2023/05/17 15:21:10 patrick Exp $
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

% set default variables
if ~exist('binlat','var') || isempty(binlat)
  binlat = 1; % latitude bin size in degrees
end
if ~exist('dawnmax','var') || isempty(dawnmax)
  dawnmax = 13;
end
if ~exist('duskmin','var') || isempty(duskmin)
  duskmin = 11;
end

% get the variables as in Li et al. 2004
mu = cosd(params.ext.oza1b);
mu0 = cosd(params.ext.sza1b);
I = params.W;

i1 = find(mu < 0);
i2 = find(mu0 < 0);
x = [params.ext.lt1b(i1);params.ext.lt1b(i2)];
xlim = [min(x),max(x)];
y = [params.ext.lat1b(i1);params.ext.lat1b(i2)];
ylim = [min(y),max(y)];

subplot(211),
scatter(params.ext.lt1b(i1),params.ext.lat1b(i1),5,params.ext.oza1b(i1))
x = [params.ext.lt1b(i1);params.ext.lt1b(i2)];
set(gca,'xlim',xlim,'ylim',ylim);
colorbar
xlabel('local time @ 1bar');
ylabel('latitutde @ 1bar');
title('Observer zenithal angle @ 1bar')

subplot(212),
scatter(params.ext.lt1b(i2),params.ext.lat1b(i2),5,params.ext.sza1b(i2))
set(gca,'xlim',xlim,'ylim',ylim);
colorbar
xlabel('local time @ 1bar');
ylabel('latitude @ 1bar');
title('Solar zenithal angle @ 1bar')

x1 = mu;
x2 = mu0;
y = I;

% conditions
lit = params.W > 0 & mu>0 & mu0 > 0; % pixel from planet surface

planet = abs(params.ext.lat1b)~=100 & ... % latitude between [-90,90] deg
         params.ext.lt1b ~= -1 & ...      % local time between [0,24] hr
         params.ext.oza1b ~= -10 & ...    % obs zenithal angle [0,90] deg
         params.ext.sza1b ~= -10;         % sol zenithal angle [0,90] deg

ltdawn = params.ext.lt1b <= dawnmax; % morning sector
ltdusk = params.ext.lt1b >= duskmin; % afternoon sector
ltovrl = params.ext.lt1b >= duskmin & params.ext.lt1b <= dawnmax; % overlap

% latitude binning
minlat = min(params.ext.lat1b(lit & planet));
maxlat = max(params.ext.lat1b(lit & planet));
nlats = fix((maxlat-minlat)/binlat);
lats = linspace(minlat, maxlat, nlats);

% calculate average value of airglow without auroral emission
glow = abs(params.ext.lat1b)<=60;
glowthrs = mean(params.W(glow))+std(params.W(glow));
fprintf(1,'Airglow threshold (mean+std)       %f\n', glowthrs);
glowthrs = median(params.W(glow))+mad(params.W(glow),1);
fprintf(1,'Airglow threshold (median+mad)     %f\n', glowthrs);
% https://uk.mathworks.com/help/matlab/ref/isoutlier.html
c = -1/(sqrt(2)*erfcinv(3/2));
glowthrs = median(params.W(glow))+c*mad(params.W(glow),1);
fprintf(1,'Airglow threshold (median+c*mad)   %f\n', glowthrs);
glowthrs = median(params.W(glow))+3*c*mad(params.W(glow),1);
fprintf(1,'Airglow threshold (median+3*c*mad) %f\n', glowthrs);
[~,~,glowthrs,~] = isoutlier(params.W(glow),'median');
fprintf(1,'Airglow threshold (isoutlier)      %f\n', glowthrs);

pause

% initialise airglow array with zeros
airglow = zeros(size(params.ext.lat1b));

% loop over latitude bins
for i = 1:nlats-1

  latmin = lats(i);
  latmax = lats(i+1);
	fprintf(1,'\nlatmin %f latmax %f\n', latmin, latmax);

  latband = planet & params.ext.lat1b <= latmax & params.ext.lat1b > latmin;
  % latitudes in morning sector
  dawn = find(lit & planet & latband & ltdawn);
  % latitudes in evening sector
  dusk = find(lit & planet & latband & ltdusk);
	% overlaping sectors
  ovrl = find(lit & planet & latband & ltovrl);

  glow = [dusk;dawn];
	if isempty(glow)
	  break; 
	end

  %glowthrs = median(params.W(glow))+mad(params.W(glow));
  glowthrs = median(params.W(glow))+c*mad(params.W(glow));
  %glowthrs = median(params.W(glow))+3*c*mad(params.W(glow));
  %[~,~,glowthrs,~] = isoutlier(params.W(glow),'median');
  fprintf(1,'Li: Airglow threshold %f\n', glowthrs);

	if exist('thrshold.mat','file')
	  glowthrs = getthrshold(.5*(lats(i)+lats(i+1)));
		fprintf(1,'Li: Airglow threshold %f (from gethrshold)\n', glowthrs);
	end

  % morning sector without aurora 
  dawnthrs = find(lit & planet & latband & ltdawn & params.W <= glowthrs);

  % evening sector without aurora
  duskthrs = find(lit & planet & latband & ltdusk & params.W <= glowthrs);

  if length(dawnthrs)>1 % morning
    % 1D fitting to x1 = mu = cosd(params.ext.oza1b)
    pdawn1 = polyfit(x1(dawnthrs),y(dawnthrs),1);
    fprintf(1,'dawn = %+f x1 %+f\n',pdawn1)
    fdawn1 = polyval(pdawn1,x1(dawn));
    % 1D fitting to x2 = mu0 = cosd(params.ext.sza1b)
    pdawn2 = polyfit(x2(dawnthrs),y(dawnthrs),1);
    fprintf(1,'dawn = %+f x2 %+f\n',pdawn2)
    fdawn2 = polyval(pdawn2,x2(dawn));
		lgdawn = {'dawn','dawn fit'};
  else
    pdawn1 = [];
    fdawn1 = NaN*ones(size(dawn));
    pdawn2 = [];
    fdawn2 = NaN*ones(size(dawn));
		lgdawn = {};
  end
  if length(duskthrs)>1 % afternoon
    % 1D fitting to x1 = mu = cosd(params.ext.oza1b)
    pdusk1 = polyfit(x1(duskthrs),y(duskthrs),1);
    fprintf(1,'dusk = %+f x1 %+f\n',pdusk1)
    fdusk1 = polyval(pdusk1,x1(dusk));
    % 1D fitting to x2 = mu0 = cosd(params.ext.sza1b)
    pdusk2 = polyfit(x2(duskthrs),y(duskthrs),1);
    fprintf(1,'dusk = %+f x2 %+f\n',pdusk2)
    fdusk2 = polyval(pdusk2,x2(dusk));
		lgdusk = {'dusk','dusk fit'};
  else
    pdusk1 = [];
    fdusk1 = NaN*ones(size(dusk));
    pdusk2 = [];
    fdusk2 = NaN*ones(size(dusk));
		lgdusk = {};
  end

  subplot(211),
  plot(x1(dawn),y(dawn),'o',x1(dawn),fdawn1,'x',...
       x1(dusk),y(dusk),'o',x1(dusk),fdusk1,'x')
  xlabel('\mu')
  ylabel('I')
	title('Li')
  % cell array {lgdawn{:},lgdusk{:}} or [lgdawn(:)',lgdusk(:)'] 
  legend([lgdawn(:)',lgdusk(:)'],'location','southeast')

  subplot(212),
  plot(x2(dawn),y(dawn),'o',x2(dawn),fdawn2,'x',...
       x2(dusk),y(dusk),'o',x2(dusk),fdusk2,'x')
  xlabel('\mu_0')
  ylabel('I')
  legend([lgdawn(:)',lgdusk(:)'],'location','southeast')

  drawnow

  if length(dawnthrs)>1 % morning
    [pdawn,gofdawn]  = fit([x1(dawnthrs),x2(dawnthrs)], y(dawnthrs),'poly11');
		disp(pdawn)
  else
    pdawn = @(x) NaN*ones(size(x,1));
  end
  if length(duskthrs)>1 % afternoon
    [pdusk,gofdusk]  = fit([x1(duskthrs),x2(duskthrs)], y(duskthrs),'poly11');
		disp(pdusk)
  else
    pdusk = @(x) NaN*ones(size(x,1));
  end

  airglow(dawn) = pdawn([x1(dawn),x2(dawn)]);
  airglow(dusk) = pdusk([x1(dusk),x2(dusk)]);
  if ~isempty(ovrl) % overlap
    x = [x1(ovrl),x2(ovrl)];
    airglow(ovrl) = 0.5*(pdawn(x)+pdusk(x));
  end
%pause

end
pause

params.li.airglow = airglow;
params.li.binlat = binlat;
params.li.dawnmax = dawnmax;
params.li.duskmin = duskmin;


clim = [min(params.W(:)),max(params.W(:))];

subplot(311)
pcolor(params.x,params.y,params.li.airglow),
shading flat
set(gca,'clim',clim);
colorbar
title(['Li ' replace(params.iFile,'_','\_')])

subplot(312)
pcolor(params.x,params.y,params.W),
shading flat
set(gca,'clim',clim);
colorbar

subplot(313)
pcolor(params.x,params.y,params.W-params.li.airglow),
shading flat
set(gca,'clim',clim);
colorbar

orient tall
print -dpdf li.pdf
pause

return
