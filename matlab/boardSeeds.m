function [S,varargout] = boardSeeds(nr,nc,ns,varargin)
% function [S[,pc]] = boardSeeds(nr,nc,ns[,'pc',pc])
% 
% string 'pc' followed by a (2 element) vector pc is optional.
% pc(1) is a percentage that indicate the size of regular tesselation
% ns = size(pc(1)*nr,pc(1)*nc) (default pc(1) = 0.1)
% pc(2) is a percentage to indicate the relative fluctuation introduced
% in the randomisation of the regular tesselation (default pc(2) = 0.075)

%
% $Id: boardSeeds.m,v 1.5 2011/03/26 17:16:55 patrick Exp $
%
% Copyright (c) 2008-2011 Patrick Guio <patrick.guio@gmail.com>
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

pc = getfield(parseArgs(struct('pc',[0.1,0.075]) , varargin{:}),'pc');

% regular tesselation with ns = pc(1)*nr x pc(1)*nc
xi = round(linspace(pc(1)/2*nc,(1-pc(1)/2)*nc, 2*round(nc*pc(1))-1));
yi = round(linspace(pc(1)/2*nr,(1-pc(1)/2)*nr, 2*round(nr*pc(1))-1));

[x, y] = meshgrid(xi,yi);
x = x(:);
y = y(:);

% initialise array S(ns,2) 
% seed s has coordinates (x,y) = S(s, 1:2) 
S = [x([1:2:end]), y([1:2:end])];
ns = size(S,1);

if pc(2), % random fluctuation of 100*pc(2) % of distance between seeds
  r = round([pc(2)*nc/4*(2*rand(ns,1)-1), pc(2)*nr/4*(2*rand(ns,1)-1)]);
  S = S + r;
end

if nargout > 1,
  varargout{1} = pc;
end
