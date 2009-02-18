function params = getDefaultparams
% function params = getDefaultparams

%
% $Id: getDefaultParams.m,v 1.1 2009/02/18 18:08:34 patrick Exp $
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
params.W = [];
params.x = [];
params.y = [];
% 
params.Wlim = [];
params.xlim = [];
params.ylim = [];

% VOISE algorithm parameters 
% Dividing
params.dividePctile = 80;
params.d2Seeds = 4;
% Merging
params.mergePctile = 60;
params.dmu = 0.2;
params.thresHoldLength = 0.3;
% Regularise
params.regMaxIter = 1;

% Diagnostics parameters
params.oDir = '../share/';
params.oMatFile = 'voise';
% Export to eps divide/merge
params.divideExport = false;
params.mergeExport = false;
% Movie parameters
params.movDiag = false;
params.movPos = [600 50 1300 1050]; % movie window size

