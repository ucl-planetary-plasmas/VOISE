function S = getCentroidSeed(VD, params, k)
% function S = getCentroidSeed(VD, params, k)

%
% $Id: getCentroidSeed.m,v 1.8 2009/08/25 15:41:59 patrick Exp $
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

% uniform weight
w = ones(size(params.W(ij)));
S = round([mean(x), mean(y)]);

% image intensity weight
if 0, % NEEDS WORK!
w = params.W(ij);
if ~any(w), % all weight equal zero
w = ones(size(params.W(ij)));
S = round([mean(x), mean(y)]);
else
S = round([sum(x.*w), sum(y.*w)]/sum(w));
end
end

if 0,% any(S<1) | any(S>size(params.W')),
fprintf(1, 'Seed %3d = (%3d, %3d), Centroid Seed = (%3d, %3d)\n', ...
	k, VD.Sx(k), VD.Sy(k), S(1), S(2));
imagesc(params.W)
axis xy
hold on
scatter(x,y,1,0.5*w);
hold off
pause
end

