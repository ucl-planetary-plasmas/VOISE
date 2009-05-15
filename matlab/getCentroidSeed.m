function S = getCentroidSeed(VD, params, k)
% function S = getCentroidSeed(VD, params, k)

%
% $Id: getCentroidSeed.m,v 1.6 2009/05/15 15:03:04 patrick Exp $
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

% find out closure of Voronoi Region associated to seed k
[ii, jj, ij] = getVRclosure(VD, k, VD.Nk{k});

x = jj;
y = ii;
w = ones(size(params.W(ij)));
S = round([mean(x), mean(y)]);
w = params.W(ij);
S = round([sum(x.*w), sum(y.*w)]/sum(w));


if 0
fprintf(1, 'Seed %3d = (%3d, %3d), Centroid Seed = (%3d, %3d)\n', ...
	k, VD.Sx(k), VD.Sy(k), S(1), S(2));
imagesc(params.x,params.y,params.W)
axis xy
hold on
scatter(x,y,5,2*w);
hold off
pause
end



