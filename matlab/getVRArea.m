function A = getVRArea(VD, sk)
% function A = getVRArea(VD, sk)

%
% $Id: getVRArea.m,v 1.3 2011/03/26 17:16:55 patrick Exp $
%
% Copyright (c) 2009-2011 Patrick Guio <patrick.guio@gmail.com>
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

% get vertices list of VR(sk)
[V,I] = getVRvertices(VD, sk);
% close the neighbour list 
is = [[1:length(I)],1];

x = [V(:,1); V(1,1)];
y = [V(:,2); V(1,2)];

if all(isfinite(x)) & all(isfinite(y))
  A = 1/2*sum(x(1:end-1).*y(2:end)-x(2:end).*y(1:end-1));
	if 0
	plot(x,y,'-x');
	for i=1:length(x)-1,
	  text(x(i), y(i), num2str(i), 'verticalalignment', 'bottom');
	end
	pause
	end
else
  A = Inf;
end
