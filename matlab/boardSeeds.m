function S = boardSeeds(nr,nc,ns,varargin)
% function S = boardSeeds(nr,nc,ns,['pc1',value],['pc2',value])
% 
% string 'pc1' followed by a value and string pc2 followed by 
% a value are optional arguments.
% pc1 is a percentage that indicate the size of regular tesselation
% ns = size(pc1*nr,pc1*nc) (default pc1 = 0.1)
% pc2 is a percentage to indicate the relative fluctuation introduced
% in the randomisation of the regular tesselation (default pc2 = 0.075)

%
% $Id: boardSeeds.m,v 1.2 2009/11/10 14:52:31 patrick Exp $
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


pc1 = getfield(parseArgs(struct('pc1',0.1)  , varargin{:}),'pc1');
pc2 = getfield(parseArgs(struct('pc2',0.075), varargin{:}),'pc2');

% regular tesselation with ns = pc1*nr x pc1*nc
xi = round(linspace(pc1/2*nc,(1-pc1/2)*nc, 2*round(nc*pc1)-1));
yi = round(linspace(pc1/2*nr,(1-pc1/2)*nr, 2*round(nr*pc1)-1));

[x, y] = meshgrid(xi,yi);
x = x(:);
y = y(:);

% initialise array S(ns,2) 
% seed s has coordinates (x,y) = S(s, 1:2) 
S = [x([1:2:end]), y([1:2:end])];
ns = size(S,1);

if pc2, % random fluctuation of 100*pc2 % of distance between seeds
  r = round([pc2*nc/4*(2*rand(ns,1)-1), pc2*nr/4*(2*rand(ns,1)-1)]);
  S = S + r;
end


