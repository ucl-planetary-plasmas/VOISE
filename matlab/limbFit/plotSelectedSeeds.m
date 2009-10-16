function plotSelectedSeeds(VD,params,fit)
% function plotSelectedSeeds(VD,params,fit)

%
% $Id: plotSelectedSeeds.m,v 1.1 2009/10/16 16:38:05 patrick Exp $
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

imagesc(params.xlim,params.ylim,params.Wo);
axis xy
axis equal
axis tight
hold on

scatter(Sx(iSelect),Sy(iSelect),LS(iSelect),LS(iSelect))
set(gca,'xlim',params.xlim,'ylim',params.ylim);

[vx,vy] = voronoi(Sx, Sy);
plot(vx,vy,'-k','LineWidth',0.5)
colorbar


ts = linspace(-180,180,50);

if length(p)==5,
  rmin = ellipse(ts,p(:).*[1;1;Rmin;Rmin;1]);
  rmax = ellipse(ts,p(:).*[1;1;Rmax;Rmax;1]);
elseif length(p)==3,
  rmin = circle(ts,p(:).*[1;1;Rmin]);
  rmax = circle(ts,p(:).*[1;1;Rmax]);
end

plot(rmin.*cosd(ts), rmin.*sind(ts), 'k-', ...
     rmax.*cosd(ts), rmax.*sind(ts), 'k-');

if length(p) == 3,

  title(sprintf('d_m=%d X_c=(%.0f,%.0f) R_0=%.0f',LSmax,p));

elseif length(p) == 5,

  title(sprintf('d_m=%d X_c=(%.0f,%.0f) a=%.0f b=%.0f tilt=%.0f',LSmax,p));

end

hold off

