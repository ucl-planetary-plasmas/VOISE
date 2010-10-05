function fit = getLimb2(VD,params,fit)
% function fit = getLimb2(VD,params,fit)

%
% $Id: getLimb2.m,v 1.4 2010/10/05 20:51:35 patrick Exp $
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

% Calculate equivalent scale length from VD as sqrt(S)
[imls, Sls] = getVDOp(VD, params.W, @(x) sqrt(length(x)));
fprintf('min(Sls) %.2f max(Sls) %.2f\n', [min(Sls), max(Sls)]);

% image axis
x = params.x;
y = params.y;

% convert from image indices to coordinates
Sx = (VD.Sx-VD.xm)/(VD.xM-VD.xm)*(max(x)-min(x))+min(x);
Sy = (VD.Sy-VD.ym)/(VD.yM-VD.ym)*(max(y)-min(y))+min(y);

% embed Sx,Sx,Sls in fit
fit.Sx  = Sx(VD.Sk);
fit.Sy  = Sy(VD.Sk);
fit.Sls = Sls;

fit = selectSeeds(fit,Sx,Sy,Sls);

% Removed unselected seeds
Sx = Sx(fit.iSelect);
Sy = Sy(fit.iSelect);
Sls = Sls(fit.iSelect);

plotSelectedSeeds(VD,params,fit);
pause


% Calculate VD statistics 
[S,xc,xy,md2s,md2c] = getVRstats(VD, params, fit.iSelect);
fprintf('min(Sls) %.2f max(Sls) %.2f\n', [min(Sls), max(Sls)]);
fprintf('min(Sls) %.2f max(Sls) %.2f\n', [min(sqrt(S)), max(sqrt(S))]);

Sls = sqrt(S);
% weight are proportional to 1/sqrt(var) 
% where var is a measure of the variance 
% or spread of the polygon
% W = sqrt(2)./Sls;
% spread around an equivalent uniform square
%W = sqrt(6)./Sls;
% spread around an equivalent uniform disc
%W = sqrt(2*pi)./Sls;
% ``exact'' estimates of the spread around the seed
W = 1./md2s;
fit = fitLimb2(fit,Sx,Sy,W);

% convert image units into pixel units
fit.pxc = (fit.p(1)-min(x))/(max(x)-min(x))*(VD.xM-VD.xm)+VD.xm;
fit.pyc = (fit.p(2)-min(y))/(max(y)-min(y))*(VD.yM-VD.ym)+VD.ym;

fprintf(1,'limb centre %8.2f, %8.2f [pixels]\n', fit.pxc,fit.pyc);

