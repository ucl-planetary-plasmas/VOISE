function thrs = getthrshold(lat)
% function thrs = getthrshold(lat)

%
% $Id: getthrshold.m,v 1.1 2021/07/09 15:22:20 patrick Exp $
%
% Copyright (c) 2021 Patrick Guio <patrick.guio@gmail.com>
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

load thrshold.mat
ii = find(LATS>-61 | LATS<-85.5);
[pt,gof]=fit(LATS(ii),THRS(ii),'smoothingspline','SmoothingParam',.1);
thrs = pt(lat);

