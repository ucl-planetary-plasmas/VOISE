function [Gx,Gy] = SobelFilter
% function [Gx,Gy] = SobelFilter
%
% see http://en.wikipedia.org/wiki/Sobel_operator
%

%
% $Id: SobelFilter.m,v 1.1 2010/11/10 15:29:13 patrick Exp $
%
% Copyright (c) 2010
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

Gx =  [-1,  0,  1; ...
       -2,  0,  2; ...
       -1,  0,  1];

Gy = -[ 1,  2,  1; ...
        0,  0,  0; ...
       -1, -2, -1];
