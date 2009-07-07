function CVD = getCentroidVD(VD, params)
% function CVD = getCentroidVD(VD, params)

%
% $Id: getCentroidVD.m,v 1.7 2009/07/07 14:16:59 patrick Exp $
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

if params.regAlgo == 2 & exist('VOISEtiming.mat','file'),
  timing = load('VOISEtiming.mat');
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
switch params.regAlgo,
  case 0, % incremental
    CVD = computeVD(nr, nc, Sc);
	case 1, % full
	  CVD = computeVDFast(nr, nc, Sc);
	case 2, % timing based
	  ns = size(Sc,1);
    tf = polyval(timing.ptVDf, ns);
	  ti = sum(polyval(timing.ptVDa,[1:ns]));
		fprintf(1,'Est. time full(%4d:%4d)/inc(%4d:%4d) %6.1f/%6.1f s\n', ...
		        1, ns, 1, ns, tf, ti);
		tStart = tic;
		if tf < ti, % full faster than incremental
		  CVD = computeVDFast(nr, nc, Sc);
		else, % incremental faster full
		  CVD = computeVD(nr, nc, Sc);
		end
		fprintf(1,'Used time %8.1f s\n', toc(tStart));
end

iter = 1;

dist = abs(CVD.Sx(CVD.Sk)-VD.Sx(Sk)) + abs(CVD.Sy(CVD.Sk)-VD.Sy(Sk));
dist2 = sqrt((CVD.Sx(CVD.Sk)-VD.Sx(Sk)).^2 + (CVD.Sy(CVD.Sk)-VD.Sy(Sk)).^2);
fprintf(1,'Iter %2d Maximum distance seed/centre-of-mass %.1f (%.1f)\n', ...
        iter, max(dist), max(dist2));

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

  switch params.regAlgo,
	  case 0, % incremental
		  CVD = computeVD(nr, nc, Sc);
		case 1, % full
		  CVD = computeVDFast(nr, nc, Sc);
		case 2, % timing based
		  ns = size(Sc,1);
	    tf = polyval(timing.ptVDf, ns);
		  ti = sum(polyval(timing.ptVDa,[1:ns]));
		  fprintf(1,'Est. time full(%4d:%4d)/inc(%4d:%4d) %6.1f/%6.1f s\n', ...
		          1, ns, 1, ns, tf, ti);
			tStart = tic;
			if tf < ti, % full faster than incremental
			  CVD = computeVDFast(nr, nc, Sc);
			else, % incremental faster full
			  CVD = computeVD(nr, nc, Sc);
			end
			fprintf(1,'Used time %8.1f s\n', toc(tStart));
	end
	iter = iter + 1;

  dist = abs(CVD.Sx(CVD.Sk)-VD.Sx(Sk)) + abs(CVD.Sy(CVD.Sk)-VD.Sy(Sk));
  dist2 = sqrt((CVD.Sx(CVD.Sk)-VD.Sx(Sk)).^2 + (CVD.Sy(CVD.Sk)-VD.Sy(Sk)).^2);
  fprintf(1,'Iter %2d Maximum distance seed/centre-of-mass %.1f (%.1f)\n', ...
	        iter, max(dist), max(dist2));

end
fprintf(1,'Centroid Voronoi Diagram computed\n')

