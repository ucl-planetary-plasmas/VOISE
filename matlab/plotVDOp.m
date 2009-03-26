function plotVDOp(VD, W, op, varargin)
% function plotVDOp(VD, W, op, varargin)

%
% $Id: plotVDOp.m,v 1.2 2009/03/26 12:00:28 patrick Exp $
%
% Copyright (c) 2008 
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

clim = [min(W(:)) max(W(:))];

%VDW = getVDOp(VD, W, @(x) median(x));
VDW = getVDOp(VD, W, op, varargin{:});
subplot(211),
imagesc(W),
axis xy,
set(gca,'clim',clim);
colorbar

subplot(212),
imagesc(VDW),
axis xy,
set(gca,'clim',clim);
colorbar

hold on
[vx,vy]=voronoi(VD.Sx(VD.Sk), VD.Sy(VD.Sk));
plot(vx,vy,'-k','LineWidth',0.5)
hold off

[op,msg] = fcnchk(op);
fprintf(1,'Image processed with %s over Voronoi Diagram\n', func2str(op))

%pause


