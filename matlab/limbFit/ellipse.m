function r=ellipse(td,p)
% function r=ellipse(td,p)
%
% Returns the polar radius of a circle wih parameters
%
% * td is an array angles [deg]
% * p is a five parameter array [xc,yc,a,b,t0]
%   where 
%     xc x-coordinate of circle center
%     yc y-coordinate of circle center
%     a  semi-major axis
%     b  semi-minor axis
%     t0 tilt angle of semi-major axis to x-axis [deg]

%
% $Id: ellipse.m,v 1.6 2023/02/01 18:39:48 patrick Exp $
%
% Copyright (c) 2022 Patrick Guio <patrick.guio@gmail.com>
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

xc = p(1); % x-coordinate of circle center
yc = p(2); % y-coordinate of circle center
a  = p(3); % semi-major axis
b  = p(4); % semi-minor axis
t0 = p(5); % tilt angle of semi-major axis to x-axis in deg

xp = a*cosd(td);
yp = b*sind(td);

Q = rot(t0);
x = zeros(size(td));
y = zeros(size(td));

for i=1:length(td),
  xy = [xc;yc]+Q*[xp(i);yp(i)];
  x(i) = xy(1);
  y(i) = xy(2);
end

r = sqrt(x.^2+y.^2);


function Q = rot(alpha)

Q = [cosd(alpha), -sind(alpha); ...
     sind(alpha), cosd(alpha)];

