function [params,IVD,DVD,MVD,CVD] = testVOISE1(nr,nc,ns,initSeeds,varargin)
% function [params,IVD,DVD,MVD,CVD] = testVOISE1(nr,nc,ns,initSeeds,varargin)
%
% example:
% [params,IVD,DVD,MVD,CVD] = testVOISE1(100,100,12,@randomSeeds)
% [params,IVD,DVD,MVD,CVD] = testVOISE1(100,100,12,@randomSeeds,...
%                            'dividePctile',90)

%
% $Id: testVOISE1.m,v 1.5 2009/11/10 17:14:05 patrick Exp $
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

global voise

% number of rows in image corresponds to y coordinate
if ~exist('nr'),
  nr = 100;
end

% number of cols in image corresponds to x coordinate
if ~exist('nc'),
  nc = 100;
end

% number of seeds 
if ~exist('ns'),
  ns = 40;
end

% load default parameters
params = getDefaultVOISEParams;

% VOISE algorithm parameters
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
params.iFile           = [voise.root '/share/testPic.mat'];
params.oDir            = [voise.root '/share/testPic/'];
params.oMatFile        = 'voise';

% diagnostics parameters
params.divideExport    = false;
params.mergeExport     = false;
params.movDiag         = false;
params.movPos          = [600 50 1300 1050]; % movie window size

% parse online parameters
params = parseArgs(params, varargin{:});

x = linspace(-2,1,nc);
y = linspace(-1,2,nr);

[X,Y] = meshgrid(x,y);

% [X(1,1)    , Y(1,1)    ] = -2  -1
% [X(1,end)  , Y(1,end)  ] =  1  -1
% [X(end,1)  , Y(end,1)  ] = -2   2
% [X(end,end), Y(end,end)] =  1   2

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
params.Wlim = [min(params.W(:)) max(params.W(:))];
params.x    = [1:size(Z,1)];
params.xlim = [min(params.x) max(params.x)];
params.y    = [1:size(Z,2)];
params.ylim = [min(params.y) max(params.y)];

if ~exist(params.oDir,'dir')
  unix(['mkdir -p ' params.oDir]);
end

% init seed of Mersenne-Twister RNG
rand('twister',10);

% run VOISE
[params,IVD,DVD,MVD,CVD] = VOISE(params, ns, initSeeds, varargin{:});

