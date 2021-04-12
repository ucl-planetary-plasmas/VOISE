function params = getHSTInfo(params)
% function params = getHSTInfo(params)

%
% $Id: getHSTInfo.m,v 1.10 2021/04/12 19:20:40 patrick Exp $
%
% Copyright (c) 2012 Patrick Guio <patrick.guio@gmail.com>
% All Rights Reserved.
%
% This program is free software; you can redistribute it and/or modify it
% under the terms of the GNU General Public License as published by the
% Free Software Foundation; either version 3 of the License, or (at your
% option) any later version.
%
% This program is distributed in the hope that it will be useful, but
% WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General
% Public License for more details.
%
% You should have received a copy of the GNU General Public License
% along with this program. If not, see <http://www.gnu.org/licenses/>.

iFile = params.iFile;

verbose = 1;

% detect whether data is from HST
TELESCOP = getFitsKeyVal(iFile,{'TELESCOP'},verbose);

if ~isempty(TELESCOP) && strcmp(TELESCOP,'HST'),
  % instrument identifier
  HST.INSTRUME = getFitsKeyVal(iFile,{'INSTRUME'},verbose);
  % proposer's target name
  HST.TARGNAME = getFitsKeyVal(iFile,{'TARGNAME'},verbose);
  % right ascension and declination of the target (deg) (J2000)
  [HST.RA_TARG,HST.DEC_TARG] = getFitsKeyVal(iFile,...
                               {'RA_TARG','DEC_TARG'},verbose);
  % UT date and time of start of first exposure
  [HST.TDATEOBS,HST.TTIMEOBS] = getFitsKeyVal(iFile,...
                                {'TDATEOBS','TTIMEOBS'},verbose);
  % exposure duration (seconds)
  HST.EXPTIME  = getFitsKeyVal(iFile,{'EXPTIME'},verbose);
  % total exposure time (seconds)
  HST.TEXPTIME = getFitsKeyVal(iFile,{'TEXPTIME'},verbose);
  if isempty(HST.TDATEOBS) & isempty(HST.TTIMEOBS),
    [HST.TDATEOBS,HST.TTIMEOBS] = getFitsKeyVal(iFile,...
                                  {'DATE-OBS','TIME-OBS'},verbose);
  end
  % start/end time (Modified Julian Time) of 1st/last exposure 
  [HST.TEXPSTRT,HST.TEXPEND] = getFitsKeyVal(iFile,...
                               {'TEXPSTRT','TEXPEND'},verbose);
  [HST.EXPSTART,HST.EXPEND] = getFitsKeyVal(iFile,...
                               {'EXPSTART','EXPEND'},verbose);
  % aperture field of view
  HST.APER_FOV = getFitsKeyVal(iFile,{'APER_FOV'},verbose);
  % plate scale (arcsec/pixel)
  if strcmp(HST.INSTRUME,'ACS'),
    HST.PLATESC = 0.025;
  else
    HST.PLATESC = getFitsKeyVal(iFile,{'PLATESC'},verbose);
  end
  % subarray axes centre point in unbinned detector pixels
  [HST.CENTERA1,HST.CENTERA2] = getFitsKeyVal(iFile,...
                                {'CENTERA1','CENTERA2'},verbose);
  % subarray axes size in unbinned detector pixels
  [HST.SIZAXIS1,HST.SIZAXIS2] = getFitsKeyVal(iFile,...
                                {'SIZAXIS1','SIZAXIS2'},verbose);
  % coordinate values at reference point
  [HST.CRVAL1,HST.CRVAL2] = getFitsKeyVal(iFile,{'CRVAL1','CRVAL2'},verbose);
  % pixel coordinates of the reference pixel
  [HST.CRPIX1,HST.CRPIX2] = getFitsKeyVal(iFile,{'CRPIX1','CRPIX2'},verbose);
  % axis type
  [HST.CTYPE1,HST.CTYPE2] = getFitsKeyVal(iFile,{'CTYPE1','CTYPE2'},verbose);
  % linear transform matrix with scale 
  % from pixel coordinate to intermediate world coordinate
  [HST.CD1_1,HST.CD1_2] = getFitsKeyVal(iFile,{'CD1_1','CD1_2'},verbose);
  [HST.CD2_1,HST.CD2_2] = getFitsKeyVal(iFile,{'CD2_1','CD2_2'},verbose);
  % inverse linear transform matrix with scale
  % from intermediate world coordinate to pixel coordinate
  if ~isempty([HST.CD1_1,HST.CD1_2;HST.CD2_1,HST.CD2_2]),
    invCD = inv([HST.CD1_1,HST.CD1_2;HST.CD2_1,HST.CD2_2]);
    HST.iCD1_1 = invCD(1,1);
    HST.iCD1_2 = invCD(1,2);
    HST.iCD2_1 = invCD(2,1);
    HST.iCD2_2 = invCD(2,2);
  end

  % ra and dec of aperture reference position
  [HST.RA_APER,HST.DEC_APER] = getFitsKeyVal(iFile,...
                               {'RA_APER','DEC_APER'},verbose);
  % offset in X/Y to subsection start 
  [HST.LTV1,HST.LTV2] = getFitsKeyVal(iFile,{'LTV1','LTV2'},verbose);
  % reciprocal of sampling rate in X/Y
  [HST.LTM1_1,HST.LTM2_2] = getFitsKeyVal(iFile,...
                            {'LTM1_1','LTM2_2'},verbose);
  % Position Angle of reference aperture center (deg)
  HST.ORIENTAT = getFitsKeyVal(iFile,{'PA_APER'},verbose);
  % position angle of image y axis (deg. e of n)
  HST.ORIENTAT = getFitsKeyVal(iFile,{'ORIENTAT'},verbose);
  % angle between sun and V1 axis (optical axis)
  HST.SUNANGLE = getFitsKeyVal(iFile,{'SUNANGLE'},verbose);

  fprintf(1,'CRPIX1 , CRPIX2   = %12.6f, %12.6f pixel\n',HST.CRPIX1,HST.CRPIX2);
  fprintf(1,'CRVAL1 , CRVAL2   = %12.6f, %12.6f deg\n',HST.CRVAL1,HST.CRVAL2);
  fprintf(1,'RA_TARG, DEC_TARG = %12.6f, %12.6f deg\n',HST.RA_TARG,HST.DEC_TARG);
  fprintf(1,'RA_APER, DEC_APER = %12.6f, %12.6f deg\n',HST.RA_APER,HST.DEC_APER);

  fprintf(1,'ORIENTAT          = %12.6f deg\n',HST.ORIENTAT);
  if ~isempty([HST.CD1_1,HST.CD1_2;HST.CD2_1,HST.CD2_2]),
    % find scaling and rotation
    CD = [HST.CD1_1,HST.CD1_2;HST.CD2_1,HST.CD2_2];
    % scaling
    s = [norm(CD(1,:));norm(CD(2,:))];
    % linear tranformation
    m = [CD(1,:)/s(1);CD(2,:)/s(2)];
    % rotation angle
    orientat = atan2(m(2,1), m(2,2))*180/pi;
    fprintf(1,'orientat CD(2,:)  = %12.6f deg\n', orientat);
	  % arcseconds per pixel
	  PIXSIZE  = s*3600;
    fprintf(1,'PIXSIZE           = %12.6f, %12.6f arcsec/pixel\n', PIXSIZE);

    % embed scaling and rotation stuff into HST 
    HST.CD = CD;
    HST.s = s;
    HST.m = m;
    HST.orientat = orientat;
    HST.PIXSIZE  = PIXSIZE;

    % inverse linear transform matrix with scale
    HST.iCD = [HST.iCD1_1,HST.iCD1_2;HST.iCD2_1,HST.iCD2_2];
  end

  % https://uk.mathworks.com/matlabcentral/answers/21209-convert-modified-julian-date
  MJD_epoch='Nov 17, 1858,00:00';
  if ~isempty(HST.EXPSTART),
	  HST.EXPSTART = datestr(HST.EXPSTART+datenum(MJD_epoch));
    fprintf(1,'EXPSTART          = %s\n', HST.EXPSTART);
	end
  if ~isempty(HST.EXPEND),
	  HST.EXPEND = datestr(HST.EXPEND+datenum(MJD_epoch));
    fprintf(1,'EXPEND            = %s\n', HST.EXPEND);
  end

  % Embed HST into params
	params.HST = HST;

end

