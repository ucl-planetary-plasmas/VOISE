function [x,y,aa,majAxis,bb,minAxis] = getEllipseAxes(a,b,c)
% function [x,y,aa,majAxis,bb,minAxis] = getEllipseAxes(a,b,c)
%
% The ellipse is assumed to be in the form
% ax^2+by^2+cxy = 1

%
% $Id: getEllipseAxes.m,v 1.3 2021/06/26 10:56:34 patrick Exp $
%
% Copyright (c) 2015 Patrick Guio <patrick.guio@gmail.com>
% All Rights Reserved.
%
% This program is free software; you can redistribute it and/or modify it
% under the terms of the GNU General Public License as published by the
% Free Software Foundation; either version 3 of the License, or (at your
% option) any later version.
%
% This program is distributed in the hope that it will be useful, but
% WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General
% Public License for more details.
%
% You should have received a copy of the GNU General Public License
% along with this program. If not, see <http://www.gnu.org/licenses/>.

% create matrix of quadratic form a x^2 + b x*y + c y^2 = 1
A = [a,c/2;c/2,b];

% get eigen vectors V and eigen values  D
[V,D] = eig(A);
%Negative V(1,1), rotates 180 degrees by changing sign
V(1,1) = - V(1,1);
%D, % eigen values are diagonal elements

% look for major axis to align with x-axis first and then rotate 
% smallest eigen value corresponds to major axis (1/sqrt(\lambda)
if D(1,1) < D(2,2), % major axis in (1,1)
  aa = 1/sqrt(D(1,1));
	tilt = atan2(V(2,1),V(1,1));
	majAxis = V(:,1);
	bb = 1/sqrt(D(2,2));
	minAxis = V(:,2);
else, % major axis in (2,2)
  aa = 1/sqrt(D(2,2));
	tilt = atan2(V(2,2),V(1,2));
	majAxis = V(:,2);
	bb = 1/sqrt(D(1,1));
	minAxis = V(:,1);
end

% check cross product of eigen vectors is +/- 1
%V(1,1)*V(2,2)-V(2,1)*V(1,2), % direct calculation
%cross([V(:,1);0],[V(:,2);0]); % 3D vector cross product 
fprintf(1,'norm V1=%.2f V2=%.2f\n', norm(majAxis), norm(minAxis));
fprintf(1,'tilt V1=%.2f V2=%.2f\n', ...
        180/pi*atan2(V(2,1),V(1,1)),180/pi*atan2(V(2,2),V(1,2)));
fprintf(1,'aa = %.2f bb = %.2f tilt = %.0f\n', aa, bb, tilt*180/pi);

t = linspace(0,2*pi,361);
if 1,
  % ellipse with major axis aligned with x
  xc = aa * cos(t);
  yc = bb * sin(t);
  % rotation by tilt
  x = xc*cos(tilt) - yc*sin(tilt);
  y = xc*sin(tilt) + yc*cos(tilt);
else
  x = aa*cos(t)*cos(tilt) - bb*sin(t)*sin(tilt);
  y = aa*cos(t)*sin(tilt) + bb*sin(t)*cos(tilt);
end

return

% visual diagnostics of the ellipse
plot(x,y,...
     majAxis(1)*[-aa;aa],majAxis(2)*[-aa;aa],...
		 minAxis(1)*[-bb;bb],minAxis(2)*[-bb;bb]);
legend('E','1','2')
axis equal
fprintf(1,'mean|a*x.^2+b*y.^2+c*x.*y-1| = %.2g\n',...
        mean(abs(a*x.^2+b*y.^2+c*x.*y-1)))
pause
