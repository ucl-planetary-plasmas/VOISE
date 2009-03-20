function [Wop, Sop] = getClosedVDOp(VD, W, op, varargin)
% function [Wop, Sop] = getClosedVDOp(VD, W, op, varargin)

%
% $Id: getClosedVDOp.m,v 1.1 2009/03/20 17:35:43 patrick Exp $
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

Wop = zeros(size(W));

if ~exist('op','var') | isempty(op), op = 'mean'; end
[op,msg] = fcnchk(op);

Swop = zeros(size(VD.Sk));
is = 1;
for s = VD.Sk', % for all seeds
  % find pixels inside the Voronoi region VR(s) or on the boundary
  ii = find(VD.Vk.lambda == s);
	% apply operator for pixels in  the Voronoi region VR(s)
	Sop(is) = op(W(ii), varargin{:});
  Wop(ii) = Sop(is);
	is = is+1;
end

