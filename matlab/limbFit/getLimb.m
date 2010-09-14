function fit = getLimb(VD,params,fit)
% function fit = getLimb(VD,params,fit)

%
% $Id: getLimb.m,v 1.2 2010/09/14 08:17:00 patrick Exp $
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

% Calculate scale length from VD 
[imls, Sls] = getVDOp(VD, params.W, @(x) sqrt(length(x)));

% image axis
x = params.x;
y = params.y;

% convert from image indices to coordinates
Sx = (VD.Sx-VD.xm)/(VD.xM-VD.xm)*(max(x)-min(x))+min(x);
Sy = (VD.Sy-VD.ym)/(VD.yM-VD.ym)*(max(y)-min(y))+min(y);

% embed Sx,Sx,Sls in fit
fit.Sx  = Sx;
fit.Sy  = Sy;
fit.Sls = Sls;

fit = selectSeeds(fit,Sx,Sy,Sls);

plotSelectedSeeds(VD,params,fit);
pause

fit = fitLimb(fit,Sx(fit.iSelect),Sy(fit.iSelect),Sls(fit.iSelect));

% convert image units into pixel units
fit.pxc = (fit.p(1)-min(x))/(max(x)-min(x))*(VD.xM-VD.xm)+VD.xm;
fit.pyc = (fit.p(2)-min(y))/(max(y)-min(y))*(VD.yM-VD.ym)+VD.ym;

fprintf(1,'limb centre %8.2f, %8.2f [pixels]\n', fit.pxc,fit.pyc);

