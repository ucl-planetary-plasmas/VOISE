function params = testReadAPISLevel2(filename)
% function params = testReadAPISLevel2(filename)
%
% Function to test APIS Level 2 images. 
%
% Examples:
%
%     * APIS HST ACS south pole image j9rlb0fxq_proc.fits
%       20/02/2007, 15:22:52-15:24:33, 100 s integration time 
%
% p = testReadAPISLevel2('../share/input/j9rlb0fxq_proc.fits');
%
%     * APIS STIS north pole image od8k1pstq_proc.fits
%       23/05/2018, 12:59:56.835, 2484.2 s integration time 
%
% p = testReadAPISLevel2('../share/input/od8k1pstq_proc.fits');
%
 

%
% $Id: testReadAPISLevel2.m,v 1.5 2023/05/18 14:49:48 patrick Exp $
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

close all

verbose = true;

params = getDefaultVOISEParams;
params.iFile = filename;
%params = getHSTInfo(params);
params = loadImage(params);

imagesc(params.x,params.y,params.W)
axis xy


if 0

info = fitsinfo(filename)
disp(info.Contents);
disp(info.PrimaryData)
disp(info.PrimaryData.Keywords);

pause

% Primary HDU
fitsdisp(filename,'mode','full','index',1)
pause

% Extended HDU
for i=2:7,  
  fitsdisp(filename,'mode','full','index',i), 
	pause, 
end

end

figure

subplot(231)
switch params.HST.APIS.HEMIS1,
  case {'South'}
  myplot(params,params.ext.lat1b,'latitude',[-90,0],-100)
	case {'North'}
  myplot(params,params.ext.lat1b,'latitude',[0,90],-100)
	otherwise
	warning('No hemisphere detected')
end

subplot(232)
myplot(params,params.ext.lt1b,'local time',[0,24],-1)

subplot(233)
myplot(params,params.ext.oza1b,'obsvr zenith angle',[0,90],-10)
[xoza,yoza] = minza(params,params.ext.oza1b);
fprintf('max oza %.1f\n',max(params.ext.oza1b(:)))

subplot(234)
myplot(params,params.ext.sza1b,'solar zenith angle',[0,90],-10)
[xsza,ysza] = minza(params,params.ext.sza1b);
fprintf('max sza %.1f\n',max(params.ext.sza1b(:)))

subplot(235)
switch params.HST.APIS.HEMIS1,
  case {'South'}
  myplot(params,params.ext.lat300km,'auroral latitude',[-90,0],-100) 
	case {'North'}
  myplot(params,params.ext.lat300km,'auroral latitude',[0,90],-100) 
	otherwise
	warning('No hemisphere detected')
end

subplot(236)
myplot(params,params.ext.lt300km,'auroral local time',[0,24],-1)

function myplot(params,z,label,clim,valNaN)

if 1
z(z==valNaN)=NaN;
[x,y]=meshgrid(params.x,params.y);
pcolor(x,y,z);
shading flat;
else
imagesc(params.x,params.y,z), 
end
axis xy, 
title(label), 
set(gca,'clim',clim),
colormap(hot)
colorbar

function [x,y] = minza(params,za)

za(za==-10)=NaN;
x = []; 
y = [];
for i=1:size(za,1),
  if ~all(~isfinite(za(i,:))) %isfinite(j)
    [mn,j] = min(za(i,:));
    x = [x,params.x(j)];
    y = [y,params.y(i)];
  end
end
hold on
plot(x,y)
hold off

