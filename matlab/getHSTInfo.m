function params = getHSTInfo(params)
% function params = getHSTInfo(params)

%
% $Id: getHSTInfo.m,v 1.13 2021/04/23 16:44:56 patrick Exp $
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

verbose = params.verbose;

HST = [];

% detect whether data is from HST
TELESCOP = getFitsKeyVal(iFile,{'TELESCOP'},verbose);

if ~isempty(TELESCOP) && strcmp(TELESCOP,'HST'), % level 1 and level 2 HST

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

  if verbose && ~isempty([HST.CRPIX1,HST.CRPIX2]),
    fprintf(1,'CRPIX1 , CRPIX2   = %12.6f, %12.6f pixel\n',...
		        HST.CRPIX1,HST.CRPIX2);
	end
  if ~isempty([HST.CRVAL1,HST.CRVAL2]),
    fprintf(1,'CRVAL1 , CRVAL2   = %12.6f, %12.6f deg\n',...
            HST.CRVAL1,HST.CRVAL2);
	end
  if ~isempty([HST.RA_TARG,HST.DEC_TARG]),
    fprintf(1,'RA_TARG, DEC_TARG = %12.6f, %12.6f deg\n',...
            HST.RA_TARG,HST.DEC_TARG);
	end
  if ~isempty([HST.RA_APER,HST.DEC_APER]),
    fprintf(1,'RA_APER, DEC_APER = %12.6f, %12.6f deg\n',...
            HST.RA_APER,HST.DEC_APER);
	end
  if ~isempty(HST.ORIENTAT),
    fprintf(1,'ORIENTAT          = %12.6f deg\n',HST.ORIENTAT);
	end
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
  HST.EPOCH_FORMAT = 'yyyy mm dd HH MM SS';
  if ~isempty(HST.EXPSTART),
	  HST.START_EPOCH = datestr(HST.EXPSTART+datenum(MJD_epoch),HST.EPOCH_FORMAT);
    fprintf(1,'START_EPOCH       = %s\n', HST.START_EPOCH);
	end
  if ~isempty(HST.EXPEND),
	  HST.END_EPOCH = datestr(HST.EXPEND+datenum(MJD_epoch),HST.EPOCH_FORMAT);
    fprintf(1,'END_EPOCH         = %s\n', HST.END_EPOCH);
  end

end
% level 2 APIS specific
SOURCE = getFitsKeyVal(iFile,{'SOURCE'},verbose);
if ~isempty(SOURCE) && strcmp(SOURCE,'APIS database'),
  fprintf(1,'Level 2 APIS fits detected.\n');

  % DATA DESCRIPTION
	HST.SOURCE = SOURCE;
  HST.APIS.TARGET1 = getFitsKeyVal(iFile,{'TARGET1'},verbose);
  HST.APIS.HEMIS1 = getFitsKeyVal(iFile,{'HEMIS1'},verbose);
  HST.APIS.HEMIS2 = getFitsKeyVal(iFile,{'HEMIS2'},verbose);

  HST.APIS.EXTEN1 = getFitsKeyVal(iFile,{'EXTEN1'},verbose);
  HST.APIS.UNIT1 = getFitsKeyVal(iFile,{'UNIT1'},verbose);

  HST.APIS.EXTEN2 = getFitsKeyVal(iFile,{'EXTEN2'},verbose);
  HST.APIS.UNIT2 = getFitsKeyVal(iFile,{'UNIT2'},verbose);

  HST.APIS.EXTEN3 = getFitsKeyVal(iFile,{'EXTEN3'},verbose);
  HST.APIS.UNIT3 = getFitsKeyVal(iFile,{'UNIT3'},verbose);

  HST.APIS.EXTEN4 = getFitsKeyVal(iFile,{'EXTEN4'},verbose);
  HST.APIS.UNIT4 = getFitsKeyVal(iFile,{'UNIT4'},verbose);

  HST.APIS.EXTEN5 = getFitsKeyVal(iFile,{'EXTEN5'},verbose);
  HST.APIS.UNIT5 = getFitsKeyVal(iFile,{'UNIT5'},verbose);

  HST.APIS.EXTEN6 = getFitsKeyVal(iFile,{'EXTEN6'},verbose);
  HST.APIS.UNIT6 = getFitsKeyVal(iFile,{'UNIT6'},verbose);

  HST.APIS.EXTEN7 = getFitsKeyVal(iFile,{'EXTEN7'},verbose);
  HST.APIS.UNIT7 = getFitsKeyVal(iFile,{'UNIT7'},verbose);

  % ORIGINAL OBSERVATION
  HST.APIS.DATEOBS = getFitsKeyVal(iFile,{'DATEOBS'},verbose);
  HST.APIS.EXP = getFitsKeyVal(iFile,{'EXP'},verbose);
  HST.APIS.PLATESC = getFitsKeyVal(iFile,{'PLATESC'},verbose);
  HST.APIS.RA_APER = getFitsKeyVal(iFile,{'RA_APER'},verbose);
  HST.APIS.DEC_APER = getFitsKeyVal(iFile,{'DEC_APER'},verbose);

  % ASTRONOMICAL EPHEMERIS AT MID-EXPOSURE
  HST.APIS.DISTE = getFitsKeyVal(iFile,{'DISTE'},verbose);
  HST.APIS.SUBELAT = getFitsKeyVal(iFile,{'SUBELAT'},verbose);
  HST.APIS.SUBELON = getFitsKeyVal(iFile,{'SUBELON'},verbose);

  HST.APIS.DISTS = getFitsKeyVal(iFile,{'DISTS'},verbose);
  HST.APIS.SUBSLAT = getFitsKeyVal(iFile,{'SUBSLAT'},verbose);
  HST.APIS.SUBSLON = getFitsKeyVal(iFile,{'SUBSLON'},verbose);

  HST.APIS.TLIGHT = getFitsKeyVal(iFile,{'TLIGHT'},verbose);
  HST.APIS.RAP = getFitsKeyVal(iFile,{'RAP'},verbose);
  HST.APIS.PHASE = getFitsKeyVal(iFile,{'PHASE'},verbose);
  HST.APIS.NP_POS = getFitsKeyVal(iFile,{'NP_POS'},verbose);
  HST.APIS.RA_TAR = getFitsKeyVal(iFile,{'RA_TAR'},verbose);
  HST.APIS.DEC_TAR = getFitsKeyVal(iFile,{'DEC_TAR'},verbose);

end

  % Embed HST into params
	params.HST = HST;


