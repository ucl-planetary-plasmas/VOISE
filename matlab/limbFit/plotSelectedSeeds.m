function plotSelectedSeeds(VD,params,fit)
% function plotSelectedSeeds(VD,params,fit)

%
% $Id: plotSelectedSeeds.m,v 1.8 2010/09/30 18:18:13 patrick Exp $
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

clf
hold on
dx = params.xlim(end)-params.xlim(1);
dy = params.ylim(end)-params.ylim(1);
dr = sqrt(dx^2+dy^2);
if ~isempty(fit.Tlim),
  for i=1:length(fit.Tlim),
	  x0 = p(1); y0 = p(2);
		x1 = x0+dr*cosd(fit.Tlim{i}(1)); y1 = y0+dr*sind(fit.Tlim{i}(1));
		x2 = x0+dr*cosd(fit.Tlim{i}(2)); y2 = y0+dr*sind(fit.Tlim{i}(2));
		fill([x0;x1;x2],[y0;y1;y2],0.7*[1,1,1]);
	end
end

scatter(Sx(iSelect),Sy(iSelect),LS(iSelect).^2,LS(iSelect),'filled')

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

  h=title(sprintf('L_M=%d C(%.1f,%.1f) R=%.1f \\epsilon(%.2f,%.2f)',...
        LSmax,p,Rmin,Rmax));

elseif length(p) == 5,

  h=title(sprintf(['L_M=%d C(%.1f,%.1f) a=%.1f b=%.1f \\alpha=%.0f ' ...
                '\\epsilon(%.2f,%.2f)'], ...
								LSmax,p,Rmin,Rmax));

end
fontsize = get(h,'FontSize');

text(0.02,0.98,sprintf('card(S)=%d',length(iSelect)),...
     'Units','Normalized',...
     'VerticalAlignment','top','HorizontalAlignment','left', ...
		 'FontSize',fontsize,'Fontweight','normal',...
		 'BackgroundColor', 0.7*[1,1,1])

if ~isempty(fit.Tlim),
  for i=1:length(fit.Tlim),
	  x0 = p(1); y0 = p(2);
		if length(p) == 3, dr = 0.8*p(3); end
		if length(p) == 5, dr = 0.4*(p(3)+p(4)); end
		x1 = x0+dr*cosd(fit.Tlim{i}(1)); y1 = y0+dr*sind(fit.Tlim{i}(1));
		text(x1,y1,sprintf('%d^\\circ',fit.Tlim{i}(1)),...
         'VerticalAlignment','bottom','HorizontalAlignment','center', ...
				 'Rotation',fit.Tlim{i}(1)-sign(fit.Tlim{i}(1))*90, ...
		     'FontSize',12,'Fontweight','normal',...
				 'BackgroundColor', 0.9*[1,1,1]);
		x2 = x0+dr*cosd(fit.Tlim{i}(2)); y2 = y0+dr*sind(fit.Tlim{i}(2));
		text(x2,y2,sprintf('%d^\\circ',fit.Tlim{i}(2)),...
         'VerticalAlignment','bottom','HorizontalAlignment','center', ...
				 'Rotation',fit.Tlim{i}(2)-sign(fit.Tlim{i}(2))*90, ...
		     'FontSize',12,'Fontweight','normal',...
				 'BackgroundColor', 0.9*[1,1,1]);
	end
end

hold off

