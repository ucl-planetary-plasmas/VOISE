function [params,IVD,DVD,MVD,CVD] = testVOISE1(varargin)
% function [params,IVD,DVD,MVD,CVD] = testVOISE1([optional args])
%
% example:
% [params,IVD,DVD,MVD,CVD] = testVOISE1(100,100,12,@randomSeeds)
% [params,IVD,DVD,MVD,CVD] = testVOISE1(100,100,12,@randomSeeds,...
%                            'dividePctile',90)

%
% $Id: testVOISE1.m,v 1.7 2009/11/12 15:22:05 patrick Exp $
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

% modify some VOISE algorithm parameters

% initialising
params.iNumSeeds       = 12;
% Dividing
params.dividePctile    = 95;
params.d2Seeds         = 4;
% Merging
params.mergePctile     = 60;
params.dmu             = 0.2;
params.thresHoldLength = 0.3;
% Regularise
params.regMaxIter      = 1;
% I/O parameters
params.iFile           = [voise.root '/share/testPic.mat'];
params.oDir            = [voise.root '/share/testPic/'];
params.oMatFile        = 'voise';
% diagnostics parameters
params.divideExport    = false;
params.mergeExport     = false;
params.movDiag         = false;
params.movPos          = [600 50 1300 1050]; % movie window size

% allow for command line parameter modifications
params = parseArgs(params, varargin{:});

try
  load(params.iFile)
  Z = double(Z);
  Z = max(Z(:))-Z;
catch
  error([params.iFile ' is not in your Matlab path\n' ...
	      'Try to run start_VOISE']);
end

% set image, axes and related
params.W    = Z;
params.x    = [1:size(Z,1)];
params.y    = [1:size(Z,2)];

params.Wlim = [min(params.W(:)) max(params.W(:))];
params.xlim = [min(params.x) max(params.x)];
params.ylim = [min(params.y) max(params.y)];

% create directory if necessary
if isunix & ~exist(params.oDir,'dir')
  unix(['mkdir -p ' params.oDir]);
end

% run VOISE
[params,IVD,DVD,MVD,CVD] = VOISE(params, varargin{:});

