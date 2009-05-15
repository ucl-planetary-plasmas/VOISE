function CVD = getCentroidVD(VD, params)
% function CVD = getCentroidVD(VD, params)

%
% $Id: getCentroidVD.m,v 1.6 2009/05/15 15:02:45 patrick Exp $
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

if params.regMaxIter < 1,
  CVD = [];
  return;
end

fprintf(1,'Computing Centroid Voronoi Diagram\n')

nr = VD.nr;
nc = VD.nc;
ns = length(VD.Sk);

% get centroid seeds
Sc = zeros(0,2);
Sk = [];
for k = VD.Sk',
  sc  = getCentroidSeed(VD, params, k);
	if isempty(find(sc(1) == Sc(:,1) & sc(2) == Sc(:,2)))
	  % in some cases two centroid seeds could be identical
		% for example here
		% o 1 1 
    % 2 x 1
		% 2 2 o
    Sc = [[Sc(:,1); sc(1)],[Sc(:,2); sc(2)]];
		Sk = [Sk; k];
	end
end
%pause
% compute centroid Voronoi Diagram 
CVD = computeVD(nr, nc, Sc);

dist = abs(CVD.Sx(CVD.Sk)-VD.Sx(Sk)) + abs(CVD.Sy(CVD.Sk)-VD.Sy(Sk));
dist2 = sqrt((CVD.Sx(CVD.Sk)-VD.Sx(Sk)).^2 + (CVD.Sy(CVD.Sk)-VD.Sy(Sk)).^2);
fprintf(1,'Maximum distance seed/centre-of-mass %.1f (%.1f)\n', ...
        max(dist), max(dist2));

iter = 1;
while max(dist) > 1e-2 & iter<params.regMaxIter, 
  % copy CVD to old VD
  VD = CVD;
  Sc = zeros(0,2);
	Sk = [];
  for k = VD.Sk',
    sc  = getCentroidSeed(VD, params, k);
		if isempty(find(sc(1) == Sc(:,1) & sc(2) == Sc(:,2)))
	    % in some cases two centroid seeds could be identical
      Sc = [[Sc(:,1); sc(1)],[Sc(:,2); sc(2)]];
			Sk = [Sk; k];
		end
  end
  %pause
  CVD = computeVD(nr, nc, Sc);
	iter = iter + 1;

  dist = abs(CVD.Sx(CVD.Sk)-VD.Sx(Sk)) + abs(CVD.Sy(CVD.Sk)-VD.Sy(Sk));
  dist2 = sqrt((CVD.Sx(CVD.Sk)-VD.Sx(Sk)).^2 + (CVD.Sy(CVD.Sk)-VD.Sy(Sk)).^2);
  fprintf(1,'Maximum distance seed/centre-of-mass %.1f (%.1f)\n', ...
	        max(dist), max(dist2));

end
fprintf(1,'Centroid Voronoi Diagram computed\n')

