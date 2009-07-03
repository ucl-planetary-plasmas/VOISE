function [Wop, Sop] = getFiniteVR(VD, W)
% function [Wop, Sop] = getFiniteVR(VD, W)

%
% $Id: getFiniteVR.m,v 1.1 2009/07/03 08:20:47 patrick Exp $
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


Wop = zeros(size(W));
Sop = zeros(size(VD.Sk));

is = 1;
for s = VD.Sk', % for all seeds
  % find pixels inside the Voronoi region VR(s) or on the boundary
  ii = find(VD.Vk.lambda == s);
	% get vertices list of VR(sk)
	[V,I] = getVRvertices(VD, s);
	if all(isfinite(V)), % finite region
	 Sop(is) = 1;
	 Wop(ii) = Sop(is);
	else
	 Sop(is) = 0;
	 Wop(ii) = Sop(is);
	end
  is = is+1;
end

