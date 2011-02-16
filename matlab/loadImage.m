function params = loadImage(params)
% function params = loadImage(params)

%
% $Id: loadImage.m,v 1.6 2011/02/16 12:44:39 patrick Exp $
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
			% overwrite pixelSize and imageOrigo dedecude from x and y
			params.pixelSize = [params.x(2)-params.x(1), params.y(2)-params.y(1)];
			params.imageOrigo = [1-params.x(1)./params.pixelSize(1),...
                           1-params.y(1)./params.pixelSize(2)];
		else
      params.x = ([1:size(params.W,2)]-params.imageOrigo(2))*params.pixelSize(2);
      params.y = ([1:size(params.W,1)]-params.imageOrigo(1))*params.pixelSize(1);
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
    params.x = ([1:size(im,2)]-params.imageOrigo(2))*params.pixelSize(2);
    params.y = ([1:size(im,1)]-params.imageOrigo(1))*params.pixelSize(1);

  else, % neither mat-file nor fits-file

    me = MException('MyFunction:fileTypeNotSupported',...
                    '%s is not a fits- nor a mat-file',params.iFile);
    throw(me);

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
