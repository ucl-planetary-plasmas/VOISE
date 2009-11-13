function [params,IVD,DVD,MVD,CVD] = VOISE(params, varargin)
% function [params,IVD,DVD,MVD,CVD] = VOISE(params, varargin)

%
% VOronoi Image SEgmentation 
%
% $Id: VOISE.m,v 1.14 2009/11/13 14:46:43 patrick Exp $
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

% miscellaneous information about VOISE
global voise

if params.logVOISE, % init diary to log VOISE run
  unix(['rm -f ' params.oDir params.oLogFile]);
  diary([params.oDir params.oLogFile])
	diary('on')
end

% start time for total time measurement of VOISE run
t = cputime;

printVOISEsetup(params);

% plot image
params = plotVOISE([], params, -1);

if params.movDiag, % init movie
  movieHandler(params, 'init');
end

[nr, nc] = size(params.W);
ns       = params.iNumSeeds;

% init seed of Mersenne-Twister RNG
rand('twister', params.RNGiseed);

if isa(params.initSeeds, 'char') | isa(params.initSeeds, 'function_handle'),
	[initSeeds, msg] = fcnchk(params.initSeeds);
  [S,params.pc] = initSeeds(nr, nc, params.iNumSeeds, varargin{:});
else
  error('initSeeds not defined or not a Function Handle');
end

if params.divideAlgo == 2 & exist([voise.root '/share/VOISEtiming.mat'],'file'),
  timing = load([voise.root '/share/VOISEtiming.mat']);
end

% save image parameters
save([params.oDir params.oMatFile], 'params'); 

% Initialise VD
fprintf(1,'*** Initialising VOISE\n')
switch params.initAlgo,
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
fprintf(1,'*** Saving VOISE results in %s\n', [params.oDir params.oMatFile]);
save([params.oDir params.oMatFile], '-append', 'CVD');
% plot
params = plotVOISE(CVD, params, 3);
if 0, 
% do not plot Voronoi diagram 
params = plotVOISE(CVD, params, 4);
end

% if movie on close movie
if params.movDiag,
  movieHandler(params, 'close');
end

t = cputime-t;
fprintf(1,'*** Total elapsed time %2d:%2d:%2d [hh:mm:ss].\n', ...
        floor(t/3600), floor(mod(t,3600)/60), floor(mod(mod(t,3660),60)));

if params.logVOISE,
	diary('off')
end
