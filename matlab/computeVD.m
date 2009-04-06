function VD = computeVD(nr, nc, S)
% function VD = computeVD(nr, nc, S)

%
% $Id: computeVD.m,v 1.2 2009/04/06 16:43:24 patrick Exp $
%
% Copyright (c) 2008 
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


VD = initVD(nr, nc, S);

ns = size(S, 1);
for k = 3:ns,
  VD = addSeedToVD(VD, S(k,:));
	if 0
  drawVD(VD);
	end
  if 0,
    fprintf(1,'.');
  end
end
if 0
fprintf(1,'\n');
end

fprintf(1,'Voronoi Diagram computed\n')
%pause
