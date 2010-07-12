function dr=dcircle2(xy,f,p,dp,func)
% function dr=dcircle2(xy,f,p,dp,func)

%
% $Id: dcircle2.m,v 1.1 2010/07/12 14:35:19 patrick Exp $
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
  fprintf(1,'calling dcircle (xc,yc,r)=(%.1f,%.1f,%.1f\n', p(1:3));
end

xc = p(1); % x-coordinate of ellipse center 
yc = p(2); % y-coordinate of ellipse center
r0 = p(3); % radius
ti = p(4:end); % angles

ni = length(xy);
ni2 = fix(ni/2);

xi = xc + xy(1:ni2);
yi = yc + xy(ni2+1:ni);

x = xc + r0*cos(ti);
y = yc + r0*sin(ti);

S = diag(sin(ti));
C = diag(cos(ti));

% dr is [drs/dxc, drs/yc, drs/dr0, 
A = [ones(size(xi)),zeros(size(xi)),cos(ti)];
B = [zeros(size(xi)),ones(size(xi)),sin(ti)];

dr = [A, -r0*S;...
      B, r0*C];


return

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


