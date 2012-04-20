function [sslat,sslong,selat,selong,sedistAU,AU2km]=computeJupiterAxis(epoch)
% function [sslat,sslong,selat,selong,sedistAU,AU2km]=computeJupiterAxis(epoch)

%
% $Id: computeJupiterAxis.m,v 1.9 2012/04/20 11:57:20 patrick Exp $
%
% Copyright (c) 2008-2012
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

% generic kernel path
persistent spiceKernelsPath

spiceKernelsPath = '/home/patrick/research/codes/spice/data/';

% Load a leapseconds kernel.
% ftp://naif.jpl.nasa.gov/pub/naif/generic_kernels/lsk/
cspice_furnsh([spiceKernelsPath 'naif0010.tls']);

% Load planetary ephemeris 
% ftp://naif.jpl.nasa.gov/pub/naif/generic_kernels/spk/planets/
cspice_furnsh([spiceKernelsPath 'de421.bsp']);

% Load satellite ephemeris 
% ftp://naif.jpl.nasa.gov/pub/naif/generic_kernels/spk/satellites/
cspice_furnsh([spiceKernelsPath 'jup267.bsp']);

% Load orientation data for planets, natural 
% satellites, the Sun, and selected asteroids
% ftp://naif.jpl.nasa.gov/pub/naif/generic_kernels/pck/
cspice_furnsh([spiceKernelsPath 'pck00010.tpc']);


% Convert string date to cspice
et = cspice_str2et(epoch);

% Get position of Sun with respect to Jupiter
target   = 'SUN';
frame    = 'IAU_JUPITER';
abcorr   = 'NONE';
observer = 'JUPITER';
% Look up the 'state' vectors and light time values 'ltime'  
% corresponding to the vector of input ephemeris time 'et'.
[state , ltime] = cspice_spkezr(target, et, frame, abcorr, observer);

% The first three entries of state contain the X, Y, Z position components.
% The final three contain the Vx, Vy, Vz velocity components.
ssposn = state(1:3);
ssdist  = norm(ssposn);
% modulo to get longitude 
sslong  = mod(atan2(ssposn(2), ssposn(1))*180/pi, 360);
sslat = 90 - acos(ssposn(3)/ssdist)*180/pi;

% and Earth
target   = 'EARTH';
frame    = 'IAU_JUPITER';
abcorr   = 'NONE';
observer = 'JUPITER';
[state , ltime] = cspice_spkezr(target, et, frame, abcorr, observer);

% Calculation of the angle between celestial north and Jupiter rotation axis 
% as seen along the line of sight Earth-Jupiter
% In IAU_JUPITER coordinates system
lineOfSight = -cspice_vhat(state(1:3));
JupiterRotAxis = [0;0;1];

% matrix to transform from Earth referential to Jupiter referential
Earth2JupiterRot = cspice_pxform('IAU_EARTH','IAU_JUPITER',et);
EarthRotAxis = cspice_vhat(Earth2JupiterRot*[0;0;1]);

% psi is the angle between the projections of the rotation axes of Earth and
% Jupiter onto the plane perpendicular to the line of sight from Earth to
% Jupiter
projJupiterRotAxis = cspice_vhat( ...
                    JupiterRotAxis-dot(JupiterRotAxis,lineOfSight)*lineOfSight);
projEarthRotAxis = cspice_vhat( ...
                    EarthRotAxis-dot(EarthRotAxis,lineOfSight)*lineOfSight);
crossProductProj = cspice_vhat(cross(projEarthRotAxis,projJupiterRotAxis));

cosa = dot(projEarthRotAxis,projJupiterRotAxis);
sina = cspice_vnorm(cross(projEarthRotAxis,projJupiterRotAxis));

if dot(crossProductProj,lineOfSight)<0,% projEarth,projJupiter anti-clockwise
  psi = acosd(cosa);
else, % clockwise
  psi = 360-acosd(cosa);
end

seposn = state(1:3);
sedist  = norm(seposn);
selong  = mod(atan2(seposn(2), seposn(1))*180/pi, 360);
selat   = 90 - acos(seposn(3)/sedist)*180/pi;


sedistAU = cspice_convrt(sedist,'KM','AU');
AU2km = cspice_convrt(1,'AU','KM');

% Transform position of Jupiter axis to Radial Tangential Normal coordinates
% Set up RTN definitions in Jupiter coordinates
% R = Sun to Jupiter unit vector
% T = (Omega x R) / | (Omega x R) | where Omega is Sun's spin axis  
% N completes the right-handed triad 

% normalised radial vector from Sun toward Jupiter 
rvec = -cspice_vhat(ssposn);

% get the matrix that transforms position vectors from Sun to Jupiter
% a specified epoch 'et'. 
% For a n-vector 'et' Sun2Jupiter is an array of dimensions (3,3,n).
Sun2Jupiter = cspice_pxform('IAU_SUN', 'IAU_JUPITER', et);

% Sun axis orientation in Sun-centred system
sunaxis = [0.0; 0.0; 1.0];
% and in frame of Jupiter
sunaxis = Sun2Jupiter * sunaxis;

% Tangential is perpendicular to radial and sun axis 
tvec = cspice_vhat(cross(sunaxis, rvec));
% Normal is perpendicular to radial and tangential
nvec = cspice_vhat(cross(rvec, tvec));

% Jupiter axis orientation in Jupiter-centred system
jupaxis = [0.0; 0.0; 1.0]; 

% in RTN
jupaxisRTN = [rvec,tvec,nvec]'*jupaxis;
if 0
jupaxis_r = sum(jupaxis.*rvec);
jupaxis_t = sum(jupaxis.*tvec);
jupaxis_n = sum(jupaxis.*nvec);
end

% get Jupiter axis / sun direction angular separation
axis_sun_ang = acos(-1.0*rvec(3))*180/pi;

if 0
fprintf(1,'%s%s\n',char(' '*ones(1,17)), ...
        'UTC   Subsol Lat  Subsol Long   Subear Lat  Subear Long');
fprintf(1,'%s %12.6f %12.6f %12.6f %12.6f\n', ...
        cspice_et2utc(et,'C',0), sslat, sslong, selat, selong);
else
fprintf(1,'Epoch                         = %s\n', cspice_et2utc(et,'C',0));
fprintf(1,'Sub-Earth latitude, longitude = %12.4f, %12.4f\n', selat, selong);
fprintf(1,'Sub-Solar latitude, longitude = %12.4f, %12.4f\n', sslat, sslong);
fprintf(1,'Earth-Jupiter distance = %.4f AU (1 AU = %.3f km)\n',sedistAU,AU2km);
end

%  It's always good form to unload kernels after use,
%  particularly in MATLAB due to data persistence.
cspice_kclear

