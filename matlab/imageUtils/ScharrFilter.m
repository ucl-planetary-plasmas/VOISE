function [Gx,Gy] = ScharrFilter
% function [Gx,Gy] = ScharrFilter
%
% see http://en.wikipedia.org/wiki/Sobel_operator
%

%
% $Id: ScharrFilter.m,v 1.2 2011/03/26 17:16:56 patrick Exp $
%
% Copyright (c) 2010-2011 Patrick Guio <patrick.guio@gmail.com>
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


Gx =  [ -3,  0,   3; ...
       -10,  0,  10; ...
        -3,  0,   3];

Gy = -[ 3,  10,   3; ...
        0,   0,   0; ...
       -3, -10,  -3];

