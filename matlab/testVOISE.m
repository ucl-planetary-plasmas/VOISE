function [params,IVD,DVD,MVD,CVD] = testVOISE(varargin)
% function [params,IVD,DVD,MVD,CVD] = testVOISE([optional args])
%
% optional arguments are any of the fieds
% examples:
% [params,IVD,DVD,MVD,CVD] = testVOISE()
% [params,IVD,DVD,MVD,CVD] = testVOISE('initSeeds',@uniformSeeds,...
%                                      'dividePctile',90)

%
% $Id: testVOISE.m,v 1.7 2009/11/12 15:16:14 patrick Exp $
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

start_VOISE

% miscellaneous information about VOISE
global voise

% load default VOISE parameters
params = getDefaultVOISEParams;

% modify some parameters

% initialising
params.iNumSeeds       = 12;
% Dividing
params.dividePctile    = 80;
params.d2Seeds         = 4;
% Merging
params.mergePctile     = 60;
params.dmu             = 0.2;
params.thresHoldLength = 0.3;
% Regularise
params.regMaxIter      = 1;
% I/O parameters
params.iFile           = [voise.root '/share/testImage.mat'];
params.oDir            = [voise.root '/share/testImage/'];
params.oMatFile        = 'voise';
% diagnostics parameters
params.divideExport    = false;
params.mergeExport     = false;
params.movDiag         = false;
params.movPos          = [600 50 1300 1050]; % movie window size

% allow command line modifications
params = parseArgs(params, varargin{:});

% load file
try
  load(params.iFile)
catch
  error([params.iFile ' is not in your Matlab path\n' ...
		     'Try to run start_VOISE']);
end

% set image, axes and related
params.W    = Z;
params.x    = x;
params.y    = y;

params.Wlim = [min(params.W(:)) max(params.W(:))];
params.xlim = [min(params.x) max(params.x)];
params.ylim = [min(params.y) max(params.y)];

% create directory if necessary
if ~exist(params.oDir,'dir')
  unix(['mkdir -p ' params.oDir]);
end


% run VOISE
[params,IVD,DVD,MVD,CVD] = VOISE(params, varargin{:});

