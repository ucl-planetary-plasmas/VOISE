function [params,VD,VD1,CVD] = VOISE(params, ns, initSeeds, varargin)
% function [params,VD,VD1,CVD] = VOISE(params, ns, initSeeds, varargin)

%
% VOronoi Image SEgmentation 
%
% $Id: VOISE.m,v 1.2 2009/02/17 14:20:12 patrick Exp $
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

%[s,w] = unix(['rm -f ' params.oDir '*.eps']);


save([params.oDir params.oMatFile], 'params'); 

if params.movDiag, % init movie
  set(gcf,'position',params.movPos);
	set(gcf,'DoubleBuffer','on');
	params.mov = avifile([params.oDir 'voise.avi'],'fps',2);
end

[nr, nc] = size(params.W);

if exist('initSeeds') & isa(initSeeds, 'function_handle'),
	[initSeeds, msg] = fcnchk(initSeeds);
  S = initSeeds(nr, nc, ns, varargin{:});
else
  error('initSeeds not defined or not a Function Handle');
end

clim = [min(params.W(:)) max(params.W(:))];
clf
subplot(111),
imagesc(params.x, params.y, params.W),
axis xy,
axis equal
%axis off
set(gca,'clim',clim);
%colorbar
set(gca,'xlim', params.xlim, 'ylim', params.ylim);
title('Original image')

exportfig(gcf,[params.oDir 'orig.eps'],'color','cmyk');

if params.movDiag,
  params.mov = addframe(params.mov, getframe(gcf,[0 0 params.movPos(3:4)]));
end

VD = computeVD(nr, nc, S);


if 0
plotVDop(VD, params.W, @(x) median(x));
%pause
end

params = plotCurrentVD(VD, params, 0);

%VD1 = [];
%CVD = [];
%return;

[VD, params] = divideVD(VD, params);
if 0
drawVD(VD);
%pause
end

save([params.oDir params.oMatFile], '-append', 'VD'); 
	
if 0
plotVDop(VD, params.W, @(x) median(x));
title('divide');
end

params = plotCurrentVD(VD, params, 1);

if 1

if ~params.movDiag, vd1 = figure; end
[VD1, params] = mergeVD(VD, params);

save([params.oDir params.oMatFile], '-append', 'VD1');

if 0
plotVDop(VD1, params.W, @(x) median(x))
title('divide+merge');
end

params = plotCurrentVD(VD1, params, 2);

if ~params.movDiag, vdc = figure; end
CVD = getCentroidVD(VD1, params);

save([params.oDir params.oMatFile], '-append', 'CVD');

if 0
plotVDop(CVD, params.W, @(x) median(x))
title('divide+merge+centroid')
end

params = plotCurrentVD(CVD, params, 3);
params = plotCurrentVD(CVD, params, 4);

if ~params.movDiag, vdl = figure; end
params = plotVDLengthScale(CVD, params);

if params.movDiag,
  params.mov = close(params.mov);
end

else, % 0

if ~params.movDiag, vdc = figure; end
CVD = getCentroidVD(VD, params);
if 0
plotVDop(CVD, params.W, @(x) median(x))
title('divide+centroid')
end

if ~params.movDiag, vd1 = figure; end
[VD1, params] = mergeVD(VVD, params);
if 0
plotVDop(VD1, params.W, @(x) median(x))
title('divide+centroid+merge');
end

end


function params = plotCurrentVD(VD, params, ic)

VDW = getVDOp(VD, params.W, @(x) median(x));

clf
subplot(111),
imagesc(params.x, params.y, VDW),
axis xy,
axis equal
%axis off
set(gca,'clim',params.Wlim);
%colorbar
%set(gca,'xlim',[VD.xm VD.xM], 'ylim', [VD.ym VD.yM]);
set(gca,'xlim', params.xlim, 'ylim', params.ylim);

if ic~=4,
hold on
[vx,vy]=voronoi(VD.Sx(VD.Sk), VD.Sy(VD.Sk));
sx = (max(params.x)-min(params.x))/(VD.xM-VD.xm);
sy = (max(params.y)-min(params.y))/(VD.yM-VD.ym);
plot((vx-VD.xm)*sx+min(params.x),(vy-VD.ym)*sy+min(params.y),'-k','LineWidth',0.5)
hold off
end

title(sprintf('card(S) = %d', length(VD.Sk)))

exportfig(gcf,[params.oDir 'phase' num2str(ic) '.eps'],'color','cmyk');

if params.movDiag,
  params.mov = addframe(params.mov, getframe(gcf,[0 0 params.movPos(3:4)]));
end

function params = plotVDLengthScale(VD, params)

VDLS = getVDOp(VD, params.W, @(x) sqrt(length(x)));
clim = [min(VDLS(VDLS>0)) max(VDLS(VDLS>0))];

clf
subplot(111),
imagesc(params.x, params.y, VDLS),
axis xy,
axis equal
%axis off
set(gca,'clim',clim);
colorbar
%set(gca,'xlim',[VD.xm VD.xM], 'ylim', [VD.ym VD.yM]);
set(gca,'xlim', params.xlim, 'ylim', params.ylim);

hold on
[vx,vy]=voronoi(VD.Sx(VD.Sk), VD.Sy(VD.Sk));
sx = (max(params.x)-min(params.x))/(VD.xM-VD.xm);
sy = (max(params.y)-min(params.y))/(VD.yM-VD.ym);
plot((vx-VD.xm)*sx+min(params.x),(vy-VD.ym)*sy+min(params.y),'-k','LineWidth',0.5)
hold off

title('Length Scale')

exportfig(gcf,[params.oDir 'ls.eps'],'color','cmyk');

if 0 & params.movPos,
  params.mov = addframe(params.mov, getframe(gcf,[0 0 params.movPos(3:4)]));
end
