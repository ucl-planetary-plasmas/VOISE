function [px, px] = getHSTradec2pixel(HST,ra,dec)
% function [px, px] = getHSTradec2pixel(HST,ra,dec)

%
% $Id: getHSTradec2pixel.m,v 1.1 2012/06/12 10:52:41 patrick Exp $
%
% Copyright (c) 2012 Patrick Guio <patrick.guio@gmail.com>
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

% reference pixel image coordinates 
rpx = HST.CRPIX1;
rpy = HST.CRPIX2;

% reference pixel ra/dec coordinates (deg)
rpra  = HST.CRVAL1;
rpdec = HST.CRVAL2;

% inverse matrix to transform from world to pixel coordinates
iCD = HST.iCD;

% planet world coordinates to pixel coordinates
px = iCD(1,1)*(ra-rpra) + iCD(1,2)*(dec-rpdec) + rpx;
yn = iCD(2,1)*(ra-rpra) + iCD(2,2)*(dec-rpdec) + rpy;

