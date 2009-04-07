function L = getVRPerimeter(VD, sk)
% function L = getVRPerimeter(VD, sk)

%
% $Id: getVRPerimeter.m,v 1.1 2009/04/07 12:27:32 patrick Exp $
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

% get vertices list of VR(sk)
[V,I] = getVRvertices(VD, sk);
% close the neighbour list 
is = [[1:length(I)],1];

x = [V(:,1); V(1,1)];
y = [V(:,2); V(1,2)];

if all(isfinite(x)) & all(isfinite(y))
  % sum(sqrt(sum((V([1:end],:)-V([2:end 1],:)).^2,2)))
  L = sum(sqrt(diff(x).^2+diff(y).^2))
	if 0
	plot(x,y,'-x');
	for i=1:length(x)-1,
	  text(x(i), y(i), num2str(i), 'verticalalignment', 'bottom');
	end
	pause
	end
else
  L = Inf;
end
