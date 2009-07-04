function params = plotDefaultVOISE(VD, params, info)
% function params = plotDefaultVOISE(VD, params, info)

%
% $Id: plotDefaultVOISE.m,v 1.1 2009/07/04 13:23:24 patrick Exp $
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

VDW = getVDOp(VD, params.W, @(x) median(x));

clf
subplot(111),
imagesc(VDW),
axis xy,
axis equal
axis off
set(gca,'clim',params.Wlim);
%colorbar
set(gca,'xlim',[VD.xm VD.xM], 'ylim', [VD.ym VD.yM]);

hold on
[vx,vy]=voronoi(VD.Sx(VD.Sk), VD.Sy(VD.Sk));
plot(vx,vy,'-k','LineWidth',0.5)
hold off

title(sprintf('card(S) = %d  (iteration %d)', length(VD.Sk), info.iter))

drawnow

if params.divideExport,
  exportfig(gcf,[params.oDir info.name num2str(info.iter) '.eps'], ...
            'color','cmyk');
end

if params.movDiag,
  params.mov = addframe(params.mov, getframe(gcf,[0 0 params.movPos(3:4)]));
end

