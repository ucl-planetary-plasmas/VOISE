function webVOISE(VOISEconf)
% function webVOISE(VOISEconf)
%
% VOISEconf is a name of a configuration file containing 
% a list of parameters for VOISE. For example
% iNumSeeds = 12
% RNGiseed = 10
% dividePctile = 80
% d2Seeds = 2
% mergePctile = 50
% dmu = 0.2
% thresHoldLength = 0.3
% regMaxIter = 2
% iFile = ../share/input/sampleint.fits
% oFile = ../share/output/sampleint/
% oMatFile = voise

% $Id: webVOISE.m,v 1.1 2010/09/03 17:20:31 patrick Exp $
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

% load default VOISE parameters
params = readVOISEconf(VOISEconf);

% load input image
params = loadImage(params);

% create directory if necessary
if ~exist(params.oDir,'dir')
    unix(['mkdir -p ' params.oDir]);
end

[params,IVD,DVD,MVD,CVD] = VOISE(params);

% some diagnostics generated in 
close all
clf
seedDist(MVD, params);
clf
plotHistHC(DVD, MVD, params);
% clf
% params = plotVDLengthScale(CVD, params);

fid = fopen([params.oDir 'CVDseeds.txt'],'w');
printSeeds(fid, CVD, params);
fclose(fid);

fid = fopen([params.oDir 'CVDneighbours.txt'],'w');
printVD(fid, CVD);
fclose(fid);

return