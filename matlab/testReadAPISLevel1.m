function testReadAPISLevel1(filename)
% function testReadAPISLevel1(filename)
%
% Function to test APIS Level 1 images.
%
% Examples:
%
% APIS HST south pole image 20/02/2007, 15:22:52-15:24:33, 100 s integration time 
% testReadAPISLevel1('../share/input/j9rlb0fxq_drz.fits'); 
%
% APIS HST north pole image 21/02/2007, 16:03:58-16:05:39, 100 s integration time 
% testReadAPISLevel1('../share/input/j9rlb1imq_drz.fits'); 
%
%

%
% $Id: testReadAPISLevel1.m,v 1.3 2023/03/28 08:34:10 patrick Exp $
%
% Copyright (c) 2021 Patrick Guio <patrick.guio@gmail.com>
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

verbose = false;
%verbose = true;

params = getDefaultVOISEParams;
params.iFile = filename;
params.verbose = verbose;
params.HSTPlanetParam = true;

params = loadImage(params);

pause

% Define right ascension-declination data set pair (in degrees)
% for the pointing direction and corresponding ephemeris times
% J2000 TDB.
raAPER = params.HST.RA_APER;
decAPER = params.HST.DEC_APER;
fprintf(1,'Aperture RA/DEC   = %f,%f deg\n', [raAPER,decAPER]);
epoch = params.HST.START_EPOCH;
loadPlanetSpiceKernels('jupiter');
et = cspice_str2et(epoch);
% Convert the RA/DEC values to radians.
d2r = cspice_rpd';
raAPER = raAPER * d2r;
decAPER = decAPER * d2r;

nels  = size(decAPER);
rad   = ones( 1,  nels(1));
z_hat = [0; 0; 1];
% Convert the angular description of the unit vectors to Cartesian.
dir   = cspice_radrec(rad, raAPER', decAPER');
% Retrieve the transformation matrix from frames J2000 to IAU_EARTH.
mat         = cspice_pxform( 'J2000', 'IAU_EARTH', et');
% Rotate the pointing direction vector into IAU_FRAME.
z = mat * dir;
[radius, raAPER, decAPER] = cspice_recrad(dir);
r2d = cspice_dpr;
raAPER  = raAPER * r2d;
decAPER = decAPER * r2d;
fprintf(1,'Aperture RA/DEC   = %f,%f deg\n', [raAPER,decAPER]);

% Calculate target RA/DEC with SPICE
target   = 'JUPITER';
frame    = 'J2000';
obsrvr   = 'EARTH';
%abcorr   = 'NONE'; % bad above
% photon leaves at et-lt and arrives at et
%abcorr   = 'LT'; % good
%abcorr   = 'LT+S'; % bad below
abcorr   = 'CN'; % good
%abcorr   = 'CN+S'; % bad below
% X* photon leaves at et and arrives at et+lt
%abcorr   = 'XLT'; % bad above
%abcorr   = 'XLT+S'; % bad above
%abcorr   = 'XCN'; % bad above
%abcorr   = 'XCN+S'; % bad above
times = [et,et+params.HST.EXPTIME];
[state,ltime] = cspice_spkezr(target,times,frame,abcorr,obsrvr);
[radJup, raJup, decJup] = cspice_recrad(state(1:3,:));
r2d = cspice_dpr;
raJup  = raJup * r2d;
decJup = decJup * r2d;
fprintf(1,'Jupiter RA/DEC 1  = %f,%f deg\n', [raJup(1),decJup(1)]);
fprintf(1,'Jupiter RA/DEC 2  = %f,%f deg\n', [raJup(2),decJup(2)]);

% get into decimal angles
hms2decdeg = @(x) (x(1)*3600+x(2)*60+x(3))/(24*3600)*360;
dms2decdeg = @(x) x(1)+x(2)/60+x(3)/3600;

[filepath,name,ext] = fileparts(filename);
switch [name,ext],
  case {'j9rlb0fxq_drz.fits'},
% from IMCCE 
%http://vo.imcce.fr/webservices/miriade/ephemph_query.php?-name=p:jupiter&-type=&-ep=%2007-02-20T15:22:52.107&-nbd=6&-step=101s&-so=1&-mime=text&-view=none&-rv=0&-anim=0&-print=1&-output=--iso,--coord:eq2000&-from=MiriadeDoc&-observer=@earth
% 20-Feb-2007 15:22:52
% RA = 17 h 2 min 46.408 s 
% DEC = -22 deg 4 min 39.53 s
% 20-Feb-2007 15:24:33
% RA = 17 h 2 min 46.446 s 
% DEC = -22 deg 4 min 39.58 s
raimcce = [hms2decdeg([17,2,46.408]), hms2decdeg([17,2,46.446])];
decimcce = [dms2decdeg(-[22,4,39.53]), dms2decdeg(-[22,4,39.58])];
  case {'j9rlb0011_drz.fits'},
% http://vo.imcce.fr/webservices/miriade/ephemph_query.php?-name=p:jupiter&-type=&-ep=%2007-02-20T15:22:52.107&-nbd=6&-step=661s&-so=1&-mime=text&-view=none&-rv=0&-anim=0&-print=1&-output=--iso,--coord:eq2000&-from=MiriadeDoc&-observer=@earth
% 20-Feb-2007 15:33:52.999
% RA = 17 h 2 min 46.655 s
% DEC = -22 deg 4 min 39.85 s
raimcce = [hms2decdeg([17,2,46.408]), hms2decdeg([17,2,46.655])];
decimcce = [dms2decdeg(-[22,4,39.53]), dms2decdeg(-[22,4,39.85])];
  case {'j9rlb1imq_drz.fits'},
% http://vo.imcce.fr/webservices/miriade/ephemph_query.php?-name=p:jupiter&-type=&-ep=%2007-02-21T16:03:58.107&-nbd=6&-step=101s&-so=1&-mime=text&-view=none&-rv=0&-anim=0&-print=1&-output=--iso,--coord:eq2000&-from=MiriadeDoc&-observer=@earth
% 21-Feb-2007 16:03:58
% RA = 17 h 3 min 19.202 s
% DEC = -22 deg 5 min 22.31 s
% 21-Feb-2007 16:05:39
% RA = 17 h 3 min 19.239 s
% DEC = -22 deg 5 min 22.35 s
raimcce = [hms2decdeg([17,3,19.202]), hms2decdeg([17,3,19.239])];
decimcce = [dms2decdeg(-[22,5,22.31]), dms2decdeg(-[22,5,22.35])];
end
if exist('raimcce','var') && exist('decimcce','var')
fprintf(1,'IMCCE RA/DEC 1    = %f,%f deg\n', raimcce(1), decimcce(1))
fprintf(1,'IMCCE RA/DEC 2    = %f,%f deg\n', raimcce(2), decimcce(2))
end

% from _proc RA_TAR, DEC_TAR
%[(raJup-ra)/params.HST.PIXSIZE(1)*3600,
%(decJup-dec)/params.HST.PIXSIZE(2)*3600]

if 0, 
[raJup,decJup]
[(raJup-raAPER)/params.HST.PIXSIZE(1)*3600,
(decJup-decAPER)/params.HST.PIXSIZE(2)*3600]
end

cspice_kclear

%imageData = fitsread(filename);%,'Info', info);
%size(imageData)
imageData = fitsread(filename,'image');%,'Info', info);
subplot(211),imagesc(imageData),axis xy
subplot(212), imagesc(params.x+1,params.y+1,params.W),axis xy
subplot(212), pcolor(params.x+1,params.y+1,params.W),shading flat

% Mesh in pixels [1,nx]x[1,ny]
[X,Y]=meshgrid(params.x+1,params.y+1);

subplot(212), pcolor(X,Y,params.W),shading flat

% Centre of image in pixels
X0 = params.HST.CRPIX1;
Y0 = params.HST.CRPIX2;
% Centre of image in AR/DEC
AR0 = params.HST.CRVAL1;
DEC0 = params.HST.CRVAL2;

% Mesh in AR/DEC
AR = params.HST.CD1_1*(X-X0)+params.HST.CD1_2*(Y-Y0)+AR0;
DEC = params.HST.CD2_1*(X-X0)+params.HST.CD2_2*(Y-Y0)+DEC0;

% check centre in pixel
xAPER = params.HST.iCD1_1*(raAPER-AR0)+params.HST.iCD1_2*(decAPER-DEC0)+X0;
yAPER = params.HST.iCD2_1*(raAPER-AR0)+params.HST.iCD2_2*(decAPER-DEC0)+Y0;

xJup = params.HST.iCD1_1*(raJup-AR0)+params.HST.iCD1_2*(decJup-DEC0)+X0;
yJup = params.HST.iCD2_1*(raJup-AR0)+params.HST.iCD2_2*(decJup-DEC0)+Y0;

if exist('raimcce','var') && exist('decimcce','var')
ximcce = params.HST.iCD1_1*(raimcce-AR0)+params.HST.iCD1_2*(decimcce-DEC0)+X0;
yimcce = params.HST.iCD2_1*(raimcce-AR0)+params.HST.iCD2_2*(decimcce-DEC0)+Y0;
end

switch [name,ext],
	case {'j9rlb0fxq_drz.fits'},
    % x=-904, y=-46 shift 
    xrenee = params.HST.CRPIX1+904;
		yrenee = params.HST.CRPIX2+46; 
  case {'j9rlb0011_drz.fits'},
    % x=-884, y=-86 shit 
    xrenee = params.HST.CRPIX1+884;
		yrenee = params.HST.CRPIX2+86;
end
if exist('xrenee','var') && exist('yrenee','var')
arRenee = params.HST.CD1_1*(xrenee-X0)+params.HST.CD1_2*(yrenee-Y0)+AR0;
decRenee = params.HST.CD2_1*(xrenee-X0)+params.HST.CD2_2*(yrenee-Y0)+DEC0;
fprintf(1,'Renee RA/DEC      = %f,%f deg\n', arRenee, decRenee)
end

subplot(212), 
clf
if 1, % in RA/DEC

h1 = plot(raJup,decJup,'ro-','MarkerSize',10);
text(raJup(1),decJup(1),'  1')
text(raJup(2),decJup(2),'  2')
Hs = h1; 
Ls = {'SPICE JUP CENTRE'};

hold on

if exist('raimcce','var') && exist('decimcce','var')
  h2 = plot(raimcce,decimcce,'go-','MarkerSize',10);
  text(raimcce(1),decimcce(1),'  1')
  text(raimcce(2),decimcce(2),'  2')
	Hs = [Hs; h2];
	Ls = [Ls(:)',{'IMCCE JUP CENTRE'}]; 
end

if exist('arRenee','var') && exist('decRenee','var')
  h3 = plot(arRenee,decRenee,'kx','MarkerSize',10);
	Hs = [Hs; h3];
	Ls = [Ls(:)',{'RENEE JUP CENTRE'}]; 
end

pcolor(AR,DEC,params.W),shading flat
h4 = plot(raAPER,decAPER,'mx','MarkerSize',10);
Hs = [Hs; h4];
Ls = [Ls(:)',{'APER CENTRE'}]; 
h5 = plot(AR0,DEC0,'mo','MarkerSize',10);
Hs = [Hs; h5];
Ls = [Ls(:)',{'IMG CENTRE'}]; 
hold off
xlabel('RA [deg]')
ylabel('DEC [deg]')
else, % in pixels
h1 = plot(xJup,yJup,'ro-','MarkerSize',10);
text(xJup(1),yJup(1),'  1')
text(xJup(2),yJup(2),'  2')
Hs = h1; 
Ls = {'SPICE JUP CENTRE'};
hold on

if exist('raimcce','var') && exist('decimcce','var')
  h2 = plot(ximcce,yimcce,'go-','MarkerSize',10);
  text(ximcce(1),yimcce(1),'  1')
  text(ximcce(2),yimcce(2),'  2')
	Hs = [Hs; h2];
	Ls = {Ls{:},'IMCCE JUP CENTRE'}; 
end

if exist('xrenee','var') && exist('yrenee','var')
  h3 = plot(xrenee,yrenee,'kx','MarkerSize',10);
	Hs = [Hs; h3];
	Ls = {Ls{:},'RENEE JUP CENTRE'}; 
end

pcolor(X,Y,params.W),shading flat
h4 = plot(xAPER,yAPER,'mx','MarkerSize',10);
Hs = [Hs; h4];
Ls = {Ls{:},'APER CENTRE'}; 
h5 = plot(X0,Y0,'mo','MarkerSize',10);
Hs = [Hs; h5];
Ls = {Ls{:},'IMG CENTRE'}; 
hold off
xlabel('pixels')
ylabel('pixels')

end % in AR/DEC/ in pixels

legend(Hs,Ls,'Location','southwest');

axis square
title(replace(filename,'_','\_'))

orient tall
print('-dpdf',replace(filename,'fits','pdf'));

if 1,
pause
else
return
end
pause off

planet = setUpSpice4Planet(params.HST);
[ss,se] = computePlanetAxis(planet.name,params.HST.START_EPOCH);

if exist('xrenee','var') && exist('yrenee','var')
pc = [xrenee,yrenee];
fprintf(1,'pc Renee = %.1f, %.1f\n', pc);
end

if exist('ximcce','var') && exist('yimcce','var')
pc = [ximcce(1),yimcce(1)];
fprintf(1,'pc IMCCE = %.1f, %.1f\n', pc);
end

pc = [xJup(1),yJup(1)];
fprintf(1,'pc Spice = %.1f, %.1f\n', pc);

epoch = params.HST.START_EPOCH;
CML = se.CML;
psi = se.psi;
orientat = params.HST.orientat;
PIXSIZE = params.HST.PIXSIZE;

plotPlanetGrid('jupiter',params,pc,epoch,CML,psi,orientat,PIXSIZE)

pause on

return

fitsdisp(filename,'mode','full','index',1);
fitsdisp(filename,'mode','full','index',2);
pause

info = fitsinfo(filename);
info
disp(info.Contents);
disp(info.PrimaryData.Keywords);

% https://uk.mathworks.com/matlabcentral/answers/21209-convert-modified-julian-date
MJD_epoch='Nov 17, 1858,00:00';

verbose=true;
dateobs = getFitsKeyVal(filename,{'DATE-OBS'},verbose)
expstart = getFitsKeyVal(filename,{'EXPSTART'},verbose)
expend = getFitsKeyVal(filename,{'EXPEND'},verbose)
exptime = getFitsKeyVal(filename,{'EXPTIME'},verbose)
datestr(expstart+datenum(MJD_epoch))
datestr(expend+datenum(MJD_epoch))

crval1 = getFitsKeyVal(filename,{'CRVAL1'},verbose)
crval2 = getFitsKeyVal(filename,{'CRVAL2'},verbose)
crpix1 = getFitsKeyVal(filename,{'CRPIX1'},verbose)
crpix2 = getFitsKeyVal(filename,{'CRPIX2'},verbose)
platesc = getFitsKeyVal(filename,{'PLATESC'},verbose)

imageData = fitsread(filename,'image','Info', info);
size(imageData)

% Fitsread cannot read PrimaryData, but can read primary
% Changed from primarydata to primary in fitsread
pData = fitsread(filename, 'primary');
size(pData)
% Fitsinfo can read PrimaryData
pDataAlt = info.PrimaryData;
size(pDataAlt)

tableData = fitsread(filename,'binarytable');
size(tableData)

