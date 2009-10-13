function dr=dcircle(t,f,p,dp,func)
% function dr=dcircle(t,f,p,dp,func)

%
% $Id: dcircle.m,v 1.1 2009/10/13 16:04:19 patrick Exp $
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

global verbose

if verbose(3),
  fprintf(1,'calling dcircle (xc,yc,r)=(%.1f,%.1f,%.1f\n', p(1:3));
end

xc = p(1); % x-coordinate of ellipse center 
yc = p(2); % y-coordinate of ellipse center
r0 = p(3); % radius

if 0,
  % circle in parametric form
  ct = cosd(t);
  st = sind(t);
	x = xc + r0*ct;
	y = yc + r0*st;

  % (r,theta) represent the circle in polar coordinates
	% t is *not* strictly equal to theta
	% therefore (r,t) does not represent the circle in polar coordinates
  % r = sqrt((xc+r0*ct).^2+(yc+r0*st).^2);
  r = f;
	theta = 180/pi*atan2(y,x);

  % dr is [dr/dxc, dr/yc, dr/dr0]
  dr = [x./r ...
        y./r ...
        (x.*ct + y.*st)./r];
else
  % sample a circle in parametric form
	ts = linspace(-180,180,7200)';
	ct = cosd(ts);
	st = sind(ts);
	xs = xc + r0*ct;
	ys = yc + r0*st;
	% and in polar coordinates
	rs = sqrt(xs.^2 + ys.^2);
	thetas = 180/pi*atan2(ys,xs);
	% drs is [drs/dxc, drs/yc, drs/dr0]
	drs = [xs./rs ...
	       ys./rs ...
				 (xs.*ct + ys.*st)./rs];

  dr = zeros(length(t),3);
	% "interpolate" dr(t) from drs(thetas)
	% interp1 cannot be used since thetas is not monotonic
	% instead find nearest approximation
	for i=1:length(t),
	  [d(i),imin] = min(abs(t(i)-thetas));
		dr(i,:) = drs(imin,:);
	end
end


end



