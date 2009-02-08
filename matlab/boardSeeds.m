function S = boardSeeds(nr,nc,ns,pc1,pc2)
% function S = boardSeeds(nr,nc,ns,pc1,pc2)

%
% $Id: boardSeeds.m,v 1.1 2009/02/08 21:07:15 patrick Exp $
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

% initialise array S(ns,2) 
% seed s has coordinates (x,y) = S(s, 1:2) 

if ~exist('pc1','var'),
	pc1 = 0.1;
end

% regular tesselation with ns = pc1*nr x pc1*nc
xi = round(linspace(pc1/2*nc,(1-pc1/2)*nc, 2*round(nc*pc1)-1));
yi = round(linspace(pc1/2*nr,(1-pc1/2)*nr, 2*round(nr*pc1)-1));


[x, y] = meshgrid(xi,yi);
x = x(:);
y = y(:);

S = [x([1:2:end]), y([1:2:end])];
ns = size(S,1);

if ~exist('pc2','var'),
	pc2 = 0.075;
end

if pc2, % random fluctuation of 100*pc2 % of distance between seeds
  r = round([pc2*nc/4*(2*rand(ns,1)-1), pc2*nr/4*(2*rand(ns,1)-1)]);
  S = S + r;
end


