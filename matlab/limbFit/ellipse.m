function r=ellipse(t,p)
% function r=ellipse(t,p)

%
% $Id: ellipse.m,v 1.3 2009/11/06 17:13:42 patrick Exp $
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

if ~isempty(verbose) & verbose(3),
  fprintf(1,'calling ellipse (xc,yc,a,b,t0)=(%.1f,%.1f,%.1f,%.1f,%.0f)\n',...
          p(1:5));
end

xc = p(1); % x-coordinate of ellipse center 
yc = p(2); % y-coordinate of ellipse center
a  = p(3); % semi-major axis
b  = p(4); % semi-minor axis
t0 = p(5); % tilt angle of semi-major axis to x-axis [deg]

% sample an ellipse in parametric form
ts = linspace(-180,180,7200)';
xs = xc + a*cosd(ts)*cosd(t0)-b*sind(ts)*sind(t0);
ys = yc + a*cosd(ts)*sind(t0)+b*sind(ts)*cosd(t0);

% and transform into polar coordinates
rs = sqrt(xs.^2 + ys.^2);
thetas = 180/pi*atan2(ys,xs);

% "interpolate" r(t) from rs(thetas)
% interp1 cannot be used since thetas is not monotonic
% instead find nearest approximation
x = zeros(size(t));
y = zeros(size(t));
r = zeros(size(t));
for i=1:length(t),
  [d(i),imin] = min(abs(t(i)-thetas));
  x(i) = xs(imin);
  y(i) = ys(imin);
  r(i) = rs(imin);
end

if 0,
  subplot(211),
  plot(x, y, 'o', r.*cosd(t), r.*sind(t), 'o')
  legend('parametric','polar');
  axis equal
  subplot(212),
  plot(t,r,'o');
end
