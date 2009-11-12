function params = getDefaultVOISEParams
% function params = getDefaultVOISEParams

%
% $Id: getDefaultVOISEParams.m,v 1.3 2009/11/12 15:13:00 patrick Exp $
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

% Image parameters
params.W    = [];
params.x    = [];
params.y    = [];
%  
params.Wlim = [];
params.xlim = [];
params.ylim = [];

% VOISE algorithm parameters 
% Initialise
params.initSeeds = @randomSeeds;
params.iNumSeeds = 12;
params.RNGiseed  = 10;
% Dividing
params.dividePctile = 80;
params.d2Seeds      = 4;
% divideAlgo choice
% 0 incremental 
% 1 full
% 2 timing based
params.divideAlgo   = 2;

% Merging
params.mergePctile     = 60;
params.dmu             = 0.2;
params.thresHoldLength = 0.3;
% mergeAlgo choice
% 0 incremental 
% 1 full
% 2 timing based
params.mergeAlgo       = 2;

% Regularise
params.regMaxIter = 1;
% regAlgo choice
% 0 incremental 
% 1 full
% 2 timing based
params.regAlgo    = 2;

% inputs/output paths
params.iFile    = [];               % input image file 
params.oDir     = [];               % output directory
params.oMatFile = 'voise.mat';      % output mat filename
params.oLogFile = 'voise.log';      % log filename
params.oMovFile = 'voise.avi';      % movie filename

% Diagnostics/Report parameters
%
% log output of VOISE run
params.logVOISE     = true;
%
% Export divide/merge/reg plot to eps
params.divideExport = false;
params.mergeExport  = false;
params.regExport    = false;
%
% Movie parameters
params.movDiag      = false;
params.movPos       = [600 50 1300 1050]; % movie window size

