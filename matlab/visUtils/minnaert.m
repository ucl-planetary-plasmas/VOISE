function params = minnaert(params,binlat,dawnmax,duskmin,np)
% function params = minnaert(params,binlat,dawnmax,duskmin,np)
%
% Compute a Minnaert-like model for the airglow as described in
% Mapping Jupiter’s Latitudinal Bands and Great Red Spot Using
% HST/WFPC2 Far-Ultraviolet Imaging, Vincent et al, Icarus, 2000
%
% params : a VOISE parameter structure with image as created by
%          functions getDefaultVOISEParams and loadImage.
% binlat : latitude bin size in degrees (default 1 deg)
% dawnmax: maximum local time for dawn sector fit (default 13h)
% duskmin: minimum local time for dusk sector fit (default 11h)
% np     : polynom order (default 3, allowed 1 and 3)


%
% $Id: minnaert.m,v 1.3 2023/05/17 16:18:14 patrick Exp $
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
if ~exist('binlat','var') || isempty(binlat),
  binlat = 1; % latitude bin size in degrees
end
if ~exist('dawnmax','var') || isempty(dawnmax),
  dawnmax = 13;
end
if ~exist('duskmin','var') || isempty(duskmin),
  duskmin = 11;
end
if ~exist('np','var') || isempty(np),
  np = 3;
end

% get the variables as in Vincent et al. 2000
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

x = log(mu.*mu0);
y = log(mu.*I);

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

switch np,
  case 3,
    fitfun = 'poly3';
    fmt = '%s = %+f x^3 %+f x^2 %+f x %+f\n';
  case 1,
    fitfun = 'poly1';
    fmt = '%s = %+f x %+f\n';
	otherwise,
    error('polynom is either of order 1 or 3')
end

LATS = 0.5*(lats(1:end-1)+lats(2:end))';
THRS = zeros(nlats-1,1);
MED = zeros(nlats-1,1);
MAD = zeros(nlats-1,1);

% loop over latitude bins
for i = 1:nlats-1,

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
	if isempty(glow), 
	  break; 
	end

  %glowthrs = median(params.W(glow))+mad(params.W(glow));
  glowthrs = median(params.W(glow))+c*mad(params.W(glow));
  %glowthrs = median(params.W(glow))+3*c*mad(params.W(glow));
  %[~,~,glowthrs,~] = isoutlier(params.W(glow),'median');
  fprintf(1,'Minnaert: Airglow threshold %f\n', glowthrs);

  MED(i) = median(params.W(glow));
  MAD(i) = mad(params.W(glow));
  THRS(i) = glowthrs;
	if exist('thrshold.mat','file'),
	  glowthrs = getthrshold(LATS(i));
    fprintf(1,'Minnaert: Airglow threshold %f (from getthrshold)\n', glowthrs);
	end

  % morning sector without aurora 
  dawnthrs = find(lit & planet & latband & ltdawn & params.W <= glowthrs);

	% evening sector without aurora
  duskthrs = find(lit & planet & latband & ltdusk & params.W <= glowthrs);

  % polyfit 
  if length(dawnthrs)>np, % morning
    pdawn = polyfit(x(dawnthrs),y(dawnthrs),np);
    fprintf(1,fmt,'dawn',pdawn)
    fdawn = polyval(pdawn,x(dawn));
    airglow(dawn) = exp(polyval(pdawn,x(dawn)))./mu(dawn);
    lgdawn = {'dawn','dawn fit'};
  else
    pdawn = [];
    fdawn = NaN*ones(size(dawn));
    lgdawn = {};
  end
  if length(duskthrs)>np, % afternoon
    pdusk = polyfit(x(duskthrs),y(duskthrs),np);
    fprintf(1,fmt,'dusk',pdusk)
    fdusk = polyval(pdusk,x(dusk));
    airglow(dusk) = exp(polyval(pdusk,x(dusk)))./mu(dusk);
    lgdusk = {'dusk','dusk fit'};
  else
    pdusk = [];
    fdusk = NaN*ones(size(dusk));
    lgdusk = {};
  end

  if 1,
  subplot(211),
	histogram(params.W(glow),'Normalization','count','BinMethod','auto')
	xlabel('I'); ylabel('count')
  else
  subplot(211)
  plot(x(dawn),y(dawn),'o',x(dawn),fdawn,'x',...
       x(dusk),y(dusk),'o',x(dusk),fdusk,'x')
  xlabel('log(\mu\mu_0)')
  ylabel('log(\mu I)')
  % cell array construction
  % by extraction {lgdaw{:},lgdusk{:}}
  % by concatenation [lgdawn(:)',lgdusk(:)']
  legend([lgdawn(:)',lgdusk(:)'],'location','southeast')
	end
	title('Minnaert')

  % using fit
  %opts = {'Robust', 'LAR'};
  %opts = {'Robust', 'Bisquare'};
  opts = {'Robust', 'Off'};
  if length(dawnthrs)>np, % morning
    [pdawn,gofdawn] = fit(x(dawnthrs), y(dawnthrs),fitfun, opts{:});
		disp(pdawn)
    airglow(dawn) = exp(pdawn(x(dawn)))./mu(dawn);
    lgdawn = {'dawn','dawn fit'};
  else
    pdawn = @(x) NaN*ones(size(x));
    lgdawn = {};
  end
  if length(duskthrs)>np, % afternoon
    [pdusk,gofdusk] = fit(x(duskthrs), y(duskthrs),fitfun, opts{:});
		disp(pdusk)
    airglow(dusk) = exp(pdusk(x(dusk)))./mu(dusk);
    lgdusk = {'dusk','dusk fit'};
  else
    pdusk = @(x) NaN*ones(size(x));
    lgdusk = {};
  end

	if ~isempty(ovrl), % overlap
	  airglow(ovrl) = 0.5*(exp(pdawn(x(ovrl)))+exp(pdusk(x(ovrl))))./mu(ovrl);
	end

  subplot(212)
  plot(x(dawn),y(dawn),'o',x(dawn),pdawn(x(dawn)),'x',...
       x(dusk),y(dusk),'o',x(dusk),pdusk(x(dusk)),'x')
  xlabel('log(\mu\mu_0)')
  ylabel('log(\mu I)')

  legend([lgdawn(:)',lgdusk(:)'],'location','southeast')
  drawnow

  %pause

end
%pause

clf
plot(LATS,MED,LATS,MAD,LATS,THRS),
xlabel('latitude bin')
ylabel('Intensity');
title('Minnaert')
legend({'MEDIAN','MAD','THRSHOLD'});

if ~exist('thrshold.mat','file'),
  fprintf(1,'*** thrshold.mat does not exists, avin!\n');
  save thrshold LATS MED MAD THRS
else
  fprintf(1,'*** thrshold.mat exists, not saving!\n');
end

pause


params.minnaert.airglow = airglow;
params.minnaert.binlat = binlat;
params.minnaert.dawnmax = dawnmax;
params.minnaert.duskmin = duskmin;
params.minnaert.np = np;
params.minnaert.fitfun = fitfun;

clim = [min(params.W(:)),max(params.W(:))];
subplot(311)
pcolor(params.x,params.y,params.minnaert.airglow), 
shading flat
set(gca,'clim',clim);
colorbar
title(['Minnaert np=' num2str(np) ' ' replace(params.iFile,'_','\_')])

subplot(312)
pcolor(params.x,params.y,params.W),
shading flat
set(gca,'clim',clim);
colorbar

subplot(313)
pcolor(params.x,params.y,params.W-params.minnaert.airglow),
shading flat
set(gca,'clim',clim);
colorbar

orient tall
print -dpdf minnaert.pdf
pause

return
