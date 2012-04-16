function [S,varargout] = poissonSeeds(nr,nc,ns,varargin)
% function [S[,pc]] = poissonSeeds(nr,nc,ns[,'pc',pc])
%
% string 'pc' followed by a value for pc is optional.
% pc is a percentage to indicate the relative fluctuation introduced
% in the randomisation of the regular tesselation (default pc = 0.02)

%
% $Id: poissonSeeds.m,v 1.6 2012/04/16 16:54:27 patrick Exp $
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

pc = getfield(parseArgs(struct('pc',0.02), varargin{:}),'pc');

% initialise array S(ns,2) 
% seed s has coordinates (x,y) = S(s, 1:2) 
% no seeds in the 100*pc % from the boundary of the image
S = round([(nc-2*pc*nc)*(rand(ns,1))) + pc*nc, ...
           (nr-2*pc*nr)*(rand(ns,1))) + pc*nr];


randraw('po',nc,[ns,1]);
randraw('po',nr,[ns,1]);

if nargout > 1,
  varargout{1} = pc;
end
