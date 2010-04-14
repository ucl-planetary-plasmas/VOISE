function varargout = printSeeds(fid, VD)
% function varargout = printSeeds(fid, VD)

%
% $Id: printSeeds.m,v 1.1 2010/04/14 09:38:37 patrick Exp $
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

s = [];

for i = VD.Sk', % for all seeds at current time
  s = [s sprintf('%3d %10d %10d', i, VD.Sx(i), VD.Sy(i))];
end

if ~isempty(fid),
	fprintf(fid, '%s', s);
end

if nargout>0,
	varargout(1) = s;
end


