function params = plotVOISE(VD, params, ic)
% function params = plotVOISE(VD, params, ic)

%
% $Id: plotVOISE.m,v 1.4 2009/11/12 16:03:17 patrick Exp $
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

clf
subplot(111),

x = params.x;
y = params.y;

if isempty(VD), % original image
	W = params.W;
else, % median operator
  W = getVDOp(VD, params.W, @(x) median(x));
end

imagesc(x, y, W),
axis xy,
axis equal
%axis off
set(gca,'clim',params.Wlim);
%colorbar
%set(gca,'xlim',[VD.xm VD.xM], 'ylim', [VD.ym VD.yM]);
set(gca,'xlim', params.xlim, 'ylim', params.ylim);

if ~isempty(VD) & ic~=4,
hold on
[vx,vy]=voronoi(VD.Sx(VD.Sk), VD.Sy(VD.Sk));
sx = (max(params.x)-min(params.x))/(VD.xM-VD.xm);
sy = (max(params.y)-min(params.y))/(VD.yM-VD.ym);
plot((vx-VD.xm)*sx+min(params.x),(vy-VD.ym)*sy+min(params.y),'-k','LineWidth',0.5)
hold off
end


if isempty(VD), % original image
  title('Original image')
  printFigure(gcf,[params.oDir 'orig.eps']);
else
  title(sprintf('card(S) = %d', length(VD.Sk)))
  printFigure(gcf,[params.oDir 'phase' num2str(ic) '.eps']);
end

if params.movDiag,
  movieHandler(params,'addframe');
end

