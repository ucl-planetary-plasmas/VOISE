function plotSelectedSeeds(VD,params,fit)
% function plotSelectedSeeds(VD,params,fit)

%
% $Id: plotSelectedSeeds.m,v 1.7 2010/09/27 17:03:44 patrick Exp $
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

p     = fit.p0;

iSelect = fit.iSelect;

Sx    = fit.Sx;
Sy    = fit.Sy;
LS    = fit.Sls;

LSmax = fit.LSmax;
Rmin  = fit.Rmin;
Rmax  = fit.Rmax;

if 0,
imagesc(params.xlim,params.ylim,params.Wo);
imagesc(params.xlim,params.ylim,params.W);
axis xy
axis equal
axis tight
hold on
end

scatter(Sx(iSelect),Sy(iSelect),LS(iSelect).^2,LS(iSelect),'filled')
hold on

[vx,vy] = voronoi(Sx, Sy);
plot(vx,vy,'-k','LineWidth',0.5)
axis equal
set(gca,'xlim',params.xlim,'ylim',params.ylim);
box on
if 1
colorbar
end


ts = linspace(-180,180,50);

if length(p)==5,
  rnom = ellipse(ts,p(:));
  rmin = ellipse(ts,p(:).*[1;1;Rmin;Rmin;1]);
  rmax = ellipse(ts,p(:).*[1;1;Rmax;Rmax;1]);
elseif length(p)==3,
  rnom = circle(ts,p(:));
  rmin = circle(ts,p(:).*[1;1;Rmin]);
  rmax = circle(ts,p(:).*[1;1;Rmax]);
end

plot(p(1), p(2), 'rx','MarkerSize',10);

h = plot(rnom.*cosd(ts), rnom.*sind(ts), 'r-', ...
     rmin.*cosd(ts), rmin.*sind(ts), 'r-', ...
     rmax.*cosd(ts), rmax.*sind(ts), 'r-');
set(h(2:3), 'LineWidth',2);

if length(p) == 3,

  h=title(sprintf('L_M=%d C(%.0f,%.0f) R=%.0f \\epsilon(%.2f,%.2f) card(S)=%d',...
        LSmax,p,Rmin,Rmax,length(iSelect)));

elseif length(p) == 5,

  h=title(sprintf(['L_M=%d C(%.0f,%.0f) a=%.0f b=%.0f \\alpha=%.0f ' ...
                '\\epsilon(%.2f,%.2f) card(S)=%d'], ...
								LSmax,p,Rmin,Rmax,length(iSelect)));

end
set(h,'FontSize',12);

hold off

