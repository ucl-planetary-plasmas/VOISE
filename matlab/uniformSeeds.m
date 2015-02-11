function [S,VDlim] = uniformSeeds(nr,nc,ns,clipping,seedfluct)
% function [S,VDlim] = uniformSeeds(nr,nc,ns,clipping,seedfluct)
% 
% clipping is defined as percentage of the image size from each edge, i.e.
% the vector of length four with [left,right,bottom,top]
% default is a 5% default clipping from all edge [left,right,bottom,top]
%
% ns = [nsx, nsy] and total is nsx * nsy
%
% seedfluct is the range of fluctuation in % of the distance between seeds
% around the nominal position

%
% $Id: uniformSeeds.m,v 1.7 2015/02/11 17:36:49 patrick Exp $
%
% Copyright (c) 2008-2012 Patrick Guio <patrick.guio@gmail.com>
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
% You should have received a copy of the GNU General Public License
% along with this program. If not, see <http://www.gnu.org/licenses/>.

if ~exist('clipping','var') || isempty(clipping),
  % 10% horizontal and 10% vertical default clipping [left,right,bottom,top]
  clipping = [5, 5, 5, 5];
end
pc = clipping/100;

if ~exist('pixfluct','var') || isempty(pixfluct),
  % 5% relative fluctuation to randomise the regular tesselation 
	pixfluct = 3;
end
pf = pixfluct/100;

% initialise array S(ns,2) 
% seed s has coordinates (x,y) = S(s, 1:2) 
% where 1 < x < nc and 1 < y < nr
% i.e. no seeds on the border of the image

xm = floor(2 + (nc-3) * pc(1));
xM = ceil(2 + (nc-3) * (1-pc(2)));
ym = floor(2 + (nr-3) * pc(3));
yM = ceil(2 + (nr-3) * (1-pc(4)));

if length(ns) == 1,
  ns = ns*ones(2,1);
end
nsx = ns(1);
nsy = ns(2);

% regular tesselation with ns = (pc(2)-pc(1))*nr x pc(1)*nc
xi = round(linspace(xm, xM, nsx));
yi = round(linspace(ym, yM, nsy));

[x, y] = meshgrid(xi,yi);

% initialise array S(ns,2) 
% seed s has coordinates (x,y) = S(s, 1:2) 
S = [x(:), y(:)];
ns = length(x(:));

if pf~=0, % random fluctuation of pixfluct % between seeds
  r = round([(xM-xm)*pf*(2*rand(ns,1)-1), (yM-ym)*pf*(2*rand(ns,1)-1)]);
  for i=1:length(r),
    s1 = S(i, :) + r(i, :);
	  if s1(1)>1 && s1(1)<nc-1 && s1(2)>1 && s1(2)<nr-1,
	    S(i, :) = s1;
		end
  end
end

% initialise VD seed limit structure
VDlim.xm = xm;
VDlim.xM = xM;
VDlim.ym = ym;
VDlim.yM = yM;

