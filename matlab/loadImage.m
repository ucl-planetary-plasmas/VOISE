function params = loadImage(params)
% function params = loadImage(params)
%
% Load an image into a VOISE parameters structure params.
%
% A default VOISE parameters structure can be generated using
% the function getDefaultVOISEParams
%
% Supported image formats supported currently include 
%
% * MATLAB mat-file format and is expected to contain
%   - at least one two-dimensional variable named 'Z' (size [nr,nc])
%     representing the image
%   - optionally two vectors named 'x' (size [1,nc]) and 'y' (size [nr,1])
%     representing the axes
%   if the axes are provided the fields 'pixelSize' and 'imageOrigo' of 
%   the VOISE parameters structure are updated accordingly, otherwise 
%   the values provided are used to generate the axes.
%
% * FITS format (Flexible Image Transport System)
%   (see http://heasarc.nasa.gov/docs/heasarc/fits.html)
%   Only one image is read. The axes are initialised using the fields
%   'pixelSize' and 'imageOrigo' of the VOISE parameters structure 
%   provided.
%
%   imageOrigo = [0, 0] means the pixel pointed by params.W(1,1) is
%   the origo

%
% $Id: loadImage.m,v 1.9 2011/02/25 14:12:28 patrick Exp $
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

try

  if strfind(params.iFile,'.mat'), % mat-file
    %   north_proj.mat is a mat file containing a polar projection of
    %   Jupiter observed by HST:
    %   Z           256x256         524288  double  image intensity
    %   x             1x256           2048  double  x-axis (# cols in Z)
    %   y           256x1             2048  double  y-axis (# rows in Z)
		if ~exist(params.iFile,'file'),
		  ME = MException('MyFunction:verifyFile', ...
                      ['Problem with file pointed to by ''params.iFile''' ...
											'\nCheck the value for params.iFile=%s' ...
                      '\nAnd/or try to run start_VOISE'],params.iFile);
      throw(ME);
		end
    im = load(params.iFile);
    % set image, axes and related
    params.W = im.Z;
		if isfield(im,'x') & isfield(im,'y'),
      params.x = im.x;
      params.y = im.y;
			% overwrite pixelSize and imageOrigo deduced from x and y
      params.pixelSize  = [diff(params.x(1:2)), diff(params.y(1:2))];
      params.imageOrigo = [-params.x(1)./params.pixelSize(1),...
                           -params.y(1)./params.pixelSize(2)];
		else
      % imageOrigo = (0,0) means params.W(1,1) is the origo
      [nr, nc] = size(params.W);
      params.x = ([0:nc-1]-params.imageOrigo(1))*params.pixelSize(1);
      params.y = ([0:nr-1]-params.imageOrigo(2))*params.pixelSize(2);
		end
		if isfield(im,'pixelUnit') & ...
       isa(im.pixelUnit,'cell') & ...
       length(im.pixelUnit) == 2,
			 params.pixelUnit = im.pixelUnit;
    end

  elseif strfind(params.iFile,'.fits'), % fits-file

		if ~exist(params.iFile,'file')
		  ME = MException('MyFunction:verifyFile', ...
                      ['Problem with file pointed to by ''params.iFile''' ...
											'\nCheck the value for params.iFile=%s' ...
                      '\nAnd/or try to run start_VOISE'],params.iFile);
      throw(ME);
		end
    im = fitsread(params.iFile);
    % set image, axes and related
    params.W = im;
    [nr, nc] = size(params.W);
    params.x = ([0:nc-1]-params.imageOrigo(1))*params.pixelSize(1);
    params.y = ([0:nr-1]-params.imageOrigo(2))*params.pixelSize(2);

  else, % neither mat-file nor fits-file

    me = MException('MyFunction:fileTypeNotSupported',...
                    '%s is not a fits- nor a mat-file',params.iFile);
    throw(me);

  end

	% ensure that image is in floating precision
	if isinteger(params.W),
	  params.W = single(params.W);
	end

  % set colour and axes limits
  if isempty(params.Wlim),
    params.Wlim = [min(params.W(:)) max(params.W(:))];
  end
  if isempty(params.xlim),
    params.xlim = [min(params.x) max(params.x)];
  end
  if isempty(params.ylim),
    params.ylim = [min(params.y) max(params.y)];
  end

catch me

  disp(['Problem when loading image file ' params.iFile '.']);
  rethrow(me);

end
