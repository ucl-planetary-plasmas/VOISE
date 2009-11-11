function [S,pc] = uniformSeeds(nr,nc,ns,varargin)
% function [S,pc] = uniformSeeds(nr,nc,ns,['pc',pc])
% 
% string 'pc' followed by a value for (two element array) pc is optional.
% pc(1) is a percentage that indicate the size of regular tesselation
% ns = size(pc(1)*nr,pc(1)*nc) (default pc(1) = 0.1)
% pc(2) is a percentage to indicate the relative fluctuation introduced
% in the randomisation of the regular tesselation (default pc(2) = 0.075)

%
% $Id: uniformSeeds.m,v 1.3 2009/11/11 17:44:02 patrick Exp $
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

pc = getfield(parseArgs(struct('pc',[0.1,0.075]) , varargin{:}),'pc');

% regular tesselation with ns = 100*pc(1)*nr x 100pc(1)*nc
xi = round(linspace(pc(1)/2*nc,(1-pc(1)/2)*nc, round(nc*pc(1))));
yi = round(linspace(pc(1)/2*nr,(1-pc(1)/2)*nr, round(nr*pc(1))));

[x, y] = meshgrid(xi,yi);

% initialise array S(ns,2) 
% seed s has coordinates (x,y) = S(s, 1:2) 
S = [x(:), y(:)];
ns = length(x(:));

if pc(2), % random fluctuation of 100*pc(2) % of distance between seeds
  r = round([pc(2)*nc/4*(2*rand(ns,1)-1), pc(2)*nr/4*(2*rand(ns,1)-1)]);
  S = S + r;
end


