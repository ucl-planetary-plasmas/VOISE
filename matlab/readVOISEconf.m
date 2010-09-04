function params = readVOISEconf(conffile)
% function params = readVOISEconf(conffile)
%
% Read a configuration file that contains
% VOISE parameters and returns a complete
% parameters structure. 
% Note that parameters not specified in the
% file are initialised to their default value.
% The syntax is 'key' = 'value' on each line.
% For example:
% iNumSeeds = 12
% RNGiseed = 10
% dividePctile = 80
% d2Seeds = 2
% mergePctile = 50
% dmu = 0.2
% thresHoldLength = 0.3
% regMaxIter = 2
% iFile = ../share/input/sampleint.fits
% oDir = ../share/output/sampleint/
% oMatFile = voise
% 
% Lines starting with # are ignored.


%
% $Id: readVOISEconf.m,v 1.2 2010/09/04 06:57:13 patrick Exp $
%
% Copyright (c) 2010
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
params = getDefaultVOISEParams();


% read VOISE configuration file
% all fields initialised by getDefaultVOISEParams 
% can be redefined in the configuration file
% The syntax is
% field = value
fid = fopen(conffile);
C = textscan(fid, '%s=%s','CommentStyle','#');
fclose(fid);
key = C{1};
val = C{2};
v = str2double(val); idx = ~isnan(v);
val(idx) = num2cell(v(idx));


for i=1:length(key),
  if isfield(params, key{i}),
    params.(key{i}) = val{i};
  end
end


