function r=circle(td,p)
% function r=circle(td,p)
%
% Returns the polar radius of a circle wih parameters
%
% * td is an array angles [deg]
% * p is a three parameter array [xc,yc,r0]
%   where 
%     xc x-coordinate of circle center
%     yc y-coordinate of circle center
%     r0 radius

%
% $Id: circle.m,v 1.6 2023/02/01 18:39:08 patrick Exp $
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
r0 = p(3); % radius

x = xc + r0*cosd(td);
y = yc + r0*sind(td);

r = sqrt(x.^2+y.^2);

