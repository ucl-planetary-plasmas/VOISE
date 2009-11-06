function dr=dellipse(t,r,p,dp,func)
% function dr=dellipse(t,r,p,dp,func)

%
% $Id: dellipse.m,v 1.4 2009/11/06 17:13:42 patrick Exp $
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
  fprintf(1,'calling dellipse (xc,yc,a,b,t0)=(%.1f,%.1f,%.1f,%.1f,%.0f)\n',...
          p(1:5));
end


xc = p(1); % x-coordinate of ellipse center 
yc = p(2); % y-coordinate of ellipse center
a  = p(3); % semi-major axis
b  = p(4); % semi-minor axis
t0 = p(5); % tilt angle of semi-major axis to x-axis [deg]

% sample an ellipse in parametric form
ts = linspace(-180,180,7200)';
ct = cosd(ts);
st = sind(ts);
xs = xc + a*ct*cosd(t0)-b*st*sind(t0);
ys = yc + a*ct*sind(t0)+b*st*cosd(t0);

% and transform into polar coordinates
rs = sqrt(xs.^2 + ys.^2);
thetas = 180/pi*atan2(ys,xs);

% drs is [drs/dxc, drs/dyc, drs/da, drs/db, drs/dt0]
drs = [xs./rs ...
       ys./rs ...
       (xs.*ct*cosd(t0) + ys.*ct*sind(t0))./rs ...
       (-xs.*st*sind(t0) + ys.*st*cosd(t0))./rs ...
       (xs.*(-a*ct*sind(t0)-b*st*cosd(t0))+ ...
        ys.*(a*ct*cosd(t0)-b*st*sind(t0)))./rs];

dr = zeros(length(t),5);
% "interpolate" dr(t) from drs(thetas)
% interp1 cannot be used since thetas is not monotonic
% instead find nearest approximation
for i=1:length(t),
  [d(i),imin] = min(abs(t(i)-thetas));
  dr(i,1:5) = drs(imin,1:5);
end


