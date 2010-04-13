function params = loadImage(params)
% function params = loadImage(params)

%
% $Id: loadImage.m,v 1.1 2010/04/13 15:39:25 patrick Exp $
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

% load file  (in this case this is a mat-file)
%   north_proj.mat is a mat file containing a polar projection of Jupiter
%   observed by HST
%   Z           256x256            524288  double  image intensity
%   x             1x256              2048  double  x-axis (nun cols in image)
%   y           256x1                2048  double  y-axis (num rows in image)
try
  im = load(params.iFile);
catch
  error([params.iFile ' is not in your Matlab path\n' ...
         'Try to run start_VOISE']);
end

% set image, axes and related
params.W = im.Z;
params.x = im.x;
params.y = im.y;

% set colour and axes limits
params.Wlim = [min(params.W(:)) max(params.W(:))];
params.xlim = [min(params.x) max(params.x)];
params.ylim = [min(params.y) max(params.y)];


