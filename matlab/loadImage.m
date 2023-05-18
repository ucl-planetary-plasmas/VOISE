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
% $Id: loadImage.m,v 1.28 2023/05/18 14:40:06 patrick Exp $
%
% Copyright (c) 2010-2012 Patrick Guio <patrick.guio@gmail.com>
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
% You should have received a copy of the GNU General Public License
% along with this program. If not, see <http://www.gnu.org/licenses/>.

try

	me = checkField(params,'iFile');
  if ~isempty(me) 
    throw(me),
  else
    iFile = params.iFile;
  end

  if contains(iFile,'.mat') % mat-file
    %   north_proj.mat is a mat file containing a polar projection of
    %   Jupiter observed by HST:
    %   Z           256x256         524288  double  image intensity
    %   x             1x256           2048  double  x-axis (# cols in Z)
    %   y           256x1             2048  double  y-axis (# rows in Z)
		me = checkiFile(params);
		if ~isempty(me), throw(me)
        end
    img = load(iFile);
    % set image, axes and related
    params.W = double(img.Z);
		if isfield(img,'x') && isfield(img,'y')
      params.x = double(img.x);
      params.y = double(img.y);
			% overwrite pixelSize and imageOrigo deduced from x and y
      params.pixelSize  = [diff(params.x(1:2)), diff(params.y(1:2))];
      params.imageOrigo = [-params.x(1)./params.pixelSize(1),...
                           -params.y(1)./params.pixelSize(2)];
		else
		  me = checkField(params,'imageOrigo'); 
		  if ~isempty(me), throw(me), end
		  me = checkField(params,'pixelSize'); 
		  if ~isempty(me), throw(me), end
      % imageOrigo = (0,0) means params.W(1,1) is the origo
      [nr, nc] = size(params.W);
      params.x = ([0:nc-1]-params.imageOrigo(1))*params.pixelSize(1);
      params.y = ([0:nr-1]-params.imageOrigo(2))*params.pixelSize(2);
    end
    if isfield(img,'pixelUnit') && ...
       isa(img.pixelUnit,'cell') && ...
       length(img.pixelUnit) == 2
      params.pixelUnit = img.pixelUnit;
    end

  elseif contains(iFile,'.fits') % fits-file

    me = checkiFile(params);
		if ~isempty(me), throw(me)
        end
		%  from [msg,msgID] = lastwarn;
    warning('off','MATLAB:imagesci:fitsinfo:unknownFormat');

    info = fitsinfo(iFile);
    opts = {'Info', info};

    % get HST fits parameters if requested
    if params.HSTFitsParam
      params = getHSTInfo(params);
    end

    ext = [];
    if ~isfield(info,'Image')
      img = squeeze(fitsread(iFile,'primary',opts{:}));
		% Information about APIS data can be found in Lamy et 2015, 
    % The Auroral Planetary Imaging and Spectroscopy (APIS) service,
    % Astronomy and Computing 11 (2015) 138–145
		% http://dx.doi.org/10.1016/j.ascom.2015.01.005
    elseif length(info.Image)==6 % APIS level 2 data _proc.fits
      % science image
      img = squeeze(fitsread(iFile,'primary',opts{:}));
      % latitude in degree at the 1-bar level
      ext.lat1b = fitsread(iFile,'image',1,opts{:});
      % local time at the 1-bar level
      ext.lt1b = fitsread(iFile,'image',2,opts{:});
      % observing zenith angle at the 1-bar level (extension 4 in Lamy 2015)
      ext.oza1b = fitsread(iFile,'image',3,opts{:});
      % solar zenith angle at the 1-bar level (extension 3 in Lamy 2015)
      ext.sza1b = fitsread(iFile,'image',4,opts{:});
      % auroral latitude at 300km above the 1-bar level
      ext.lat300km = fitsread(iFile,'image',5,opts{:});
      % auroral local time at 300km above the 1-bar level
      ext.lt300km = fitsread(iFile,'image',6,opts{:});
    elseif length(info.Image)==3 % APIS level 1 data
      % original calibrated image (correction for dark current,
      % flat-fielding, correction for geometric distortion)
      % *_drz.fits => calibrated ACS,  units = electrons s−1 pix−1,
      % The (dark) line of bad pixels visible was replaced by values
      % linearly interpolated from neighbor pixels.
      % *_x2d.fits => calibrated STIS, units = counts pix−1
      img = squeeze(fitsread(iFile,'image',1,opts{:})); 
      % weight(?) image
      ext.wht = squeeze(fitsread(iFile,'image',2,opts{:})); 
      % context(?) image
      ext.ctx = squeeze(fitsread(iFile,'image',3,opts{:})); 
    end
    warning('on','MATLAB:imagesci:fitsinfo:unknownFormat');

    % set image, axes and related
    params.W = img;
    params.ext = ext;

    if isfield(params,'HST') && ...
       ~isempty(params.HST) && ...
       params.HSTPlanetParam
      me = checkSpice; 
      if ~isempty(me), throw(me), end
      params = getHSTPlanetParams(params);
    end

		me = checkField(params,'imageOrigo'); 
		if ~isempty(me), throw(me), end

		me = checkField(params,'pixelSize'); 
		if ~isempty(me), throw(me), end

    [nr, nc] = size(params.W);
    params.x = ([0:nc-1]-params.imageOrigo(1))*params.pixelSize(1);
    params.y = ([0:nr-1]-params.imageOrigo(2))*params.pixelSize(2);

if 0
    % pixel coordinates (indices j)
    [Xj,Yj] = meshgrid(1:nc, 1:nr);

    HST = params.HST;
    % reference pixel image coordinates 
    rpx = HST.CRPIX1;
    rpy = HST.CRPIX2;
    % reference pixel ra/dec coordinates (deg)
    rpra  = HST.CRVAL1;
    rpdec = HST.CRVAL2;

    % pixel coordinates to world coordinates (ra/dec in deg) (indices i)
    [Xi,Yi] = getHSTpixel2radec(HST,Xj,Yj);

    % planet world coordinates (ra/dec) to pixel coordinates
    planet = params.planet;
    [planet.pxc,planet.pyc] = getHSTradec2pixel(HST,planet.ra,planet.dec);
    params.planet = planet;

if 0
    figure
    pcolor(Xj,Yj,log10(abs(params.W))); shading flat; 
    hold on
    plot(rpx,rpy,'ko','markersize',5)
    plot(planet.pxc,planet.pyc,'kx','markersize',5)
    hold off

    figure
    if 1
      % convert to arcsec and set relative ra/dec to ref pix
      [Xi,Yi] = getHSTabs2relRadec(HST,Xi,Yi);
      [rpra,rpdec] = getHSTabs2relRadec(HST,rpra,rpdec);
      [planet.ra,planet.dec] = getHSTabs2relRadec(HST,planet.ra,planet.dec);
      xlbl = 'ra/ref. pixel [arcsec]';
      ylbl = 'dec/ref.pixel [arcsec]';
    else
      % absolute ra/dec in deg
      xlbl = 'ra [deg]';
      ylbl = 'dec [deg]';
    end
    pcolor(Xi,Yi,log10(abs(params.W))); shading flat; 
    hold on
    plot(rpra,rpdec,'ko','markersize',5)
    plot(planet.ra,planet.dec,'kx','markersize',5)
    hold off
    xlabel(xlbl);
    ylabel(ylbl);
end
end

  else % neither mat-file nor fits-file

    if ~isempty(iFile)
      me = MException('MyFunction:fileTypeNotSupported',...
                      '%s is not a fits- nor a mat-file',iFile);
    else
		  me = MException('MyFunction:fileEmpty',...
                      'field iFile is empty');
      throw(me);
    end

  end

	% ensure that image is in floating precision
  if isinteger(params.W)
    params.W = single(params.W);
  end

  % apply misc. filtering if required
	params = preprocessImage(params);

  % set colour and axes limits
  if isempty(params.Wlim)
    params.Wlim = [min(params.W(:)) max(params.W(:))];
  end
  if isempty(params.xlim)
    params.xlim = [min(params.x) max(params.x)];
  end
  if isempty(params.ylim)
    params.ylim = [min(params.y) max(params.y)];
  end

catch me

  disp('Problem when loading image file.');
  rethrow(me);

end


function string = mydeblank(string)

string = deblank(fliplr(deblank(fliplr(string))));

function me = checkSpice()

me = [];
if isempty(which('mice')) || exist('mice') ~= 3
  me = MException('MyFunction:verifySpice', ...
                  ['Problem with Mice/Spice installation' ...
                  '\nCheck your installation']);
end

function me = checkiFile(params)

me = [];
if ~exist(params.iFile,'file')
  me = MException('MyFunction:verifyFile', ...
                  ['Problem with file pointed to by ''params.iFile''' ...
                  '\nCheck the value for params.iFile=%s' ...
                  '\nAnd/or try to run start_VOISE'],params.iFile);
end

function me = checkField(params,field)

me = [];
if ~isfield(params,field)
  names = fieldnames(params);
	strnames = '';
	if ~isempty(names)
	  strnames = '(';
    for i=1:length(names)-1 
      strnames = [strnames names{i} ', ']; 
      return
    end
	  strnames = [strnames names{end} ')'];
	end
  me = MException('MyFunction:verifyParams', ...
                  ['Problem with ''params'' structure fields ''%s''' ...
									'\nNon-existing fields ''%s''' ...
                  '\nYou should run getDefaultVOISEParams'], strnames, field);
end
