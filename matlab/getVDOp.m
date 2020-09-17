function [Wop, Sop] = getVDOp(VD, W, op, varargin)
% function [Wop, Sop] = getVDOp(VD, W, op, varargin)

%
% $Id: getVDOp.m,v 1.9 2020/09/17 06:30:24 patrick Exp $
%
% Copyright (c) 2008-2012 Patrick Guio <patrick.guio@gmail.com>
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

Wop = zeros(size(W));

if ~exist('op','var') || isempty(op), op = 'mean'; end

switch op,
  case 'median',
    metricID = 1; 
		op = @(x) median(x);
  case 'mean',
	  metricID = 2; 
		op = @(x) mean(x);
  case 'range',
	  metricID = 3; 
		op = @(x) (max(x)-min(x));
  case 'sqrtLen',
	  metricID = 4; 
		op = @(x) sqrt(length(x));
  case 'homogeneousCriteria', 
	  metricID = 5; 
		op = @(x,y) (max(x)-min(x))/y;
  case 'stdDev', 
	  metricID = 6; 
		op = @(x,y) y*std(x); 
	otherwise, error('op not recognized');
end

[op,msg] = fcnchk(op);

if exist('getVDOpBatch')==3, % W has to be double!
	[Wop, Sop] = getVDOpBatch(VD, W, metricID, varargin{:});
else
  Sop = zeros(size(VD.Sk));
  is = 1;
  for s = VD.Sk', % for all seeds
    % find pixels inside the Voronoi region VR(s)
    ii = find(VD.Vk.lambda == s & VD.Vk.v == 0);
	  % apply operator for pixels in  the Voronoi region VR(s)
	  Sop(is) = op(W(ii), varargin{:});
    Wop(ii) = Sop(is);
	  is = is+1;
  end
end

% find pixels not in any Voronoi region but on boundaries
ii = find(VD.Vk.v == 1);

if 0,
  % rescale the value for these points 
  Wop(ii) = min(Sop(:))-(max(Sop(:))-min(Sop(:)));
else
  Wop(ii) = NaN;
end

