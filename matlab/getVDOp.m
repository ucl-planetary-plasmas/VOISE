function [Wop, Sop] = getVDOp(VD, W, op, varargin)
% function [Wop, Sop] = getVDOp(VD, W, op, varargin)

%
% $Id: getVDOp.m,v 1.2 2009/03/18 15:49:58 patrick Exp $
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
  % find pixels in the Voronoi region VR(s)
  ii = find(VD.Vk.lambda == s & VD.Vk.v == 0);
	% apply operator for pixels in  the Voronoi region VR(s)
	Sop(is) = op(W(ii), varargin{:});
  Wop(ii) = Sop(is);
	is = is+1;
end

% find pixels not in any Voronoi region 
ii = find(VD.Vk.v == 1);
if 0,
% rescale the value for these points 
Wop(ii) = min(Sop(:))-(max(Sop(:))-min(Sop(:)));
else
Wop(ii) = NaN;
end
  
