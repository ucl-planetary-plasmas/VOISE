function testReadAPISLevel2(filename)
% function testReadAPISLevel2(filename)
%
% Function to test APIS Level 1 images. 
%
% Example:
%
% APIS HST south pole image 20/02/2007, 15:22:52-15:24:33, 100 s integration time 
% testReadAPISproc('../share/input/j9rlb0fxq_proc.fits')

%
% $Id: testReadAPISLevel2.m,v 1.2 2023/03/28 08:34:10 patrick Exp $
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

info = fitsinfo(filename);

if 0
info
disp(info.Contents);
disp(info.PrimaryData)
disp(info.PrimaryData.Keywords);

pause

% Primary HDU
fitsdisp(filename,'mode','full','index',1)
pause

% Extended HDU
for i=2:7,  fitsdisp(filename,'mode','full','index',i), pause, end

end

figure

subplot(231)
plot(params,params.ext.lat1b,'latitude',[-90,0])

subplot(232)
plot(params,params.ext.lt1b,'local time',[0,24])

subplot(233)
plot(params,params.ext.oza1b,'obsvr zenithal angle',[0,90])

subplot(234)
plot(params,params.ext.sza1b,'solar zenithal angle',[0,90])

subplot(235)
plot(params,params.ext.lat300km,'auroral latitude',[-90,0]) 

subplot(236)
plot(params,params.ext.lt300km,'auroral local time',[0,24])

function plot(params,z,label,clim)
imagesc(params.x,params.y,z), 
axis xy, 
title(label), 
set(gca,'clim',clim),
colorbar
