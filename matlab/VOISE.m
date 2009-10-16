function [params,IVD,DVD,MVD,CVD] = VOISE(params, ns, initSeeds, varargin)
% function [params,IVD,DVD,MVD,CVD] = VOISE(params, ns, initSeeds, varargin)

%
% VOronoi Image SEgmentation 
%
% $Id: VOISE.m,v 1.6 2009/10/16 13:46:32 patrick Exp $
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

%[s,w] = unix(['rm -f ' params.oDir '*.eps']);

% save image parameters
save([params.oDir params.oMatFile], 'params'); 
% plot image
params = plotVOISE([], params, -1);

if params.movDiag, % init movie
  set(gcf,'position',params.movPos);
	set(gcf,'DoubleBuffer','on');
	params.mov = avifile([params.oDir 'voise.avi'],'fps',2);
end

[nr, nc] = size(params.W);

if exist('initSeeds') & isa(initSeeds, 'function_handle'),
	[initSeeds, msg] = fcnchk(initSeeds);
  S = initSeeds(nr, nc, ns, varargin{:});
else
  error('initSeeds not defined or not a Function Handle');
end

if params.divideAlgo == 2 & exist('VOISEtiming.mat','file'),
  timing = load('VOISEtiming.mat');
end

% Initialise VD
fprintf(1,'*** Initialising VOISE\n')
switch params.divideAlgo,
  case 0, % incremental
    IVD = computeVD(nr, nc, S);
	case 1, % full
	  IVD = computeVDFast(nr, nc, S);
	case 2, % timing based
	  tf = polyval(timing.ptVDf, ns);
		ti = sum(polyval(timing.ptVDa, [1:ns]));
		fprintf(1,'Est. time full(%4d:%4d)/inc(%4d:%4d) %6.1f/%6.1f s ', ...
		        1, ns, 1, ns, tf, ti);
		tStart = tic;
		if tf < ti, % full faster than incremental
		  IVD = computeVDFast(nr, nc, S);
		else, % incremental faster than full
		  IVD = computeVD(nr, nc, S);
		end
		fprintf(1,'(Used %6.1f s)\n', toc(tStart));
end
fprintf(1,'*** Initialising completed.\n')

% save 
save([params.oDir params.oMatFile], '-append', 'IVD'); 
% plot 
params = plotVOISE(IVD, params, 0);

% Dividing phase 
[DVD, params] = divideVD(IVD, params);
% save 
save([params.oDir params.oMatFile], '-append', 'DVD'); 
% plot
params = plotVOISE(DVD, params, 1);

% if movie on do not change figure
if ~params.movDiag, vd1 = figure; end

% Merging phase
[MVD, params] = mergeVD(DVD, params);
% save 
save([params.oDir params.oMatFile], '-append', 'MVD');
% plot
params = plotVOISE(MVD, params, 2);

% if movie on do not change figure
if ~params.movDiag, vdc = figure; end

% Regularisation phase
CVD = getCentroidVD(MVD, params);
% save
save([params.oDir params.oMatFile], '-append', 'CVD');
% plot
params = plotVOISE(CVD, params, 3);
if 0, 
% do not plot Voronoi diagram 
params = plotVOISE(CVD, params, 4);
end

% if movie on close movie
if params.movDiag,
  params.mov = close(params.mov);
end

