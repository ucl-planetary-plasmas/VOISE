function [params,IVD,DVD,MVD,CVD] = testVOISE1(nr,nc,ns,initSeeds,varargin)
% function [params,IVD,DVD,MVD,CVD] = testVOISE1(nr,nc,ns,initSeeds,varargin)
%
% example:
% [params,IVD,DVD,MVD,CVD] = testVOISE1(100,100,12,@randomSeeds)

%
% $Id: testVOISE1.m,v 1.4 2009/11/06 15:43:15 patrick Exp $
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

x = linspace(-2,1,nc);
y = linspace(-1,2,nr);

[X,Y] = meshgrid(x,y);

% [X(1,1)    , Y(1,1)    ] = -2  -1
% [X(1,end)  , Y(1,end)  ] =  1  -1
% [X(end,1)  , Y(end,1)  ] = -2   2
% [X(end,end), Y(end,end)] =  1   2

load ../share/testPic
Z = double(Z);
Z = max(Z(:))-Z;

% load default parameters
params = getDefaultVOISEParams;

% set image, axes and related
params.W = Z;
params.Wlim = [min(params.W(:)) max(params.W(:))];
params.x = [1:size(Z,1)];
params.xlim = [min(params.x) max(params.x)];
params.y = [1:size(Z,2)];
params.ylim = [min(params.y) max(params.y)];

params.dividePctile = 80;
params.d2Seeds = 4;
params.mergePctile = 60;
params.dmu = 0.2;
params.thresHoldLength = 0.3;
params.regMaxIter = 1;

params.oDir = '../share/';
params.oMatFile = 'voise';
params.divideExport = false;
params.mergeExport = false;
params.movDiag = false;
params.movPos=[600 50 1300 1050]; % movie window size

% init seed of Mersenne-Twister RNG
rand('twister',10);

close all

t = cputime;
[params,IVD,DVD,MVD,CVD] = VOISE(params, ns, initSeeds, varargin{:});
t = cputime-t;
fprintf(1,'Elapsed time %2d:%2d:%2d\n', ...
  floor(t/3600), floor(mod(t,3600)/60), floor(mod(mod(t,3660),60)));

