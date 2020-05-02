function SaturnRingDemo
% function SaturnRingDemo
%
% Demo of planetary rings visualisation with a Saturn HST/ACS image
% j9rls2011_drz_F125LP.fits from 2007 JAN 13 06:38:11

%
% $Id: SaturnRingDemo.m,v 1.1 2020/05/02 17:30:33 patrick Exp $
%
% Copyright (c) 2008-2020 Patrick Guio <patrick.guio@gmail.com>
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

% load VOISE parameters
params = readVOISEconf('../share/SaturnRingDemo.dat');

% load image
params = loadImage(params);

% print VOISE set up
% printVOISEsetup(params);

% HST keywords and parameters
HST = params.HST; 

pc       = [HST.CRPIX1, HST.CRPIX2];        % planet centre from HST (pixels)
pc       = [724, 591];                      % planet centre from eyes 
epoch    = [HST.TDATEOBS ' ' HST.TTIMEOBS]; % time of observation
orientat = HST.ORIENTAT;                    % orientation angle 
PIXSIZE  = HST.s*3600;                      % arcseconds per pixel
CML      = [];
psi      = [];

% common arguments
args = {'saturn', params, pc, epoch, CML, psi, orientat, PIXSIZE};

% rings specifications
[ringsNames,ringsSpecs] = getPlanetRings('saturn');

close all
pause('off');

% plot planet grid without rings
figure, plotPlanetGrid(args{:},'no rings')

% plot planet grid with automatic rings
figure, plotPlanetGrid(args{:})

% plot planet grid with specific rings only
figure, plotPlanetGrid(args{:},ringsSpecs(1,:))

% plot only rings (no grid) with automatic rings
figure, plotPlanetGrid(args{:},'no grid')

% plot only rings (no grid) with specific radii
figure, plotPlanetGrid(args{:},'no grid',ringsSpecs(1:2,:))

