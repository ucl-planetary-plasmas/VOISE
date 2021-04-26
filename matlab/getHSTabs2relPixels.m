function [x,y] = getHSTabs2relPixels(HST,x,y)
% function [x,y] = getHSTabs2relPixels(HST,x,y)
%
% Transform absolute x,y pixel coordinates from [1,nx],[1,ny] into relative
% pixel aith respect to reference pixel in arcsec
% HST is a structure containing necessary HST fits header information 
% provided by function getHSTInfo()

%
% $Id: getHSTabs2relPixels.m,v 1.1 2021/04/26 13:07:18 patrick Exp $
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

% reference pixel ra/dec coordinates (deg)
rpx = HST.CRPIX1;
rpy = HST.CRPIX2;

% scale factor from deg to arcsec
s   = HST.PIXSIZE;  % arcsec/pixel

if length(s) == 2,
  x = s(1)*(x - rpx);
  y = s(2)*(y - rpy);
elseif length(s) == 1,
  x = s*(x - rpx);
  y = s*(y - rpy);
end

