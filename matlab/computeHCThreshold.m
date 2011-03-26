function [WD,SD,WHC,SHC,HCThreshold] = computeHCThreshold(VD, params, Pctile)
% function [WD,SD,WHC,SHC,HCThreshold] = computeHCThreshold(VD, params, Pctile)

%
% $Id: computeHCThreshold.m,v 1.2 2011/03/26 17:16:55 patrick Exp $
%
% Copyright (c) 2008-2011 Patrick Guio <patrick.guio@gmail.com>
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

% D(sk) = max(x)-min(x), x in VR(sk)
[WD, SD] = getVDOp(VD, params.W, @(x) max(x)-min(x));
maxD = max(SD);

% Homogeneity function
[WHC, SHC] = getVDOp(VD, params.W, ...
                     @(x) homogeneousCriteria(x, max(SD)));

% Homogeneity dynamic threshold
HCThreshold = prctile(SHC, Pctile);

fprintf(1,'Homogeneity Criteria min=%.2f max=%.2f Q(%.0f)=%.2f\n', ...
        [min(SHC), max(SHC)], Pctile, HCThreshold)
