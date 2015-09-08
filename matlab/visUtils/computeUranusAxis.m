function [sslat,sslong,selat,selong,CML,psi,sedistAU,AU2km]=computeUranusAxis(epoch)
% function [sslat,sslong,selat,selong,CML,psi,sedistAU,AU2km]=computeUranusAxis(epoch)

%
% $Id: computeUranusAxis.m,v 1.5 2015/09/08 14:51:41 patrick Exp $
%
% Copyright (c) 20012
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

loadPlanetSpiceKernels('uranus');

% Convert string date to cspice
et = cspice_str2et(epoch);

% Get position of Sun with respect to Uranus
target   = 'SUN';
frame    = 'IAU_URANUS';
abcorr   = 'NONE';
observer = 'URANUS';
% Look up the 'state' vectors and light time values 'ltime'  
% corresponding to the vector of input ephemeris time 'et'.
[state , ltime] = cspice_spkezr(target, et, frame, abcorr, observer);

% The first three entries of state contain the X, Y, Z position components.
% The final three contain the Vx, Vy, Vz velocity components.
ssposn = state(1:3);
ssdist  = norm(ssposn);
% modulo to get longitude 
%sslong  = mod(atan2(ssposn(2), ssposn(1))*180/pi, 360);
sslong  = atan2(ssposn(2), ssposn(1))*cspice_dpr;
sslat = 90 - acos(ssposn(3)/ssdist)*cspice_dpr;
%ssdist,[sslong,sslat]

% converts rectangular coordinates to latitudinal coordinates.
[radius, lon, lat] = cspice_reclat(ssposn);
%radius,[lon lat]*cspice_dpr

% converts rectangular coordinates to planetographic coordinates.
radii = cspice_bodvrd('URANUS','RADII',3);
re = radii(1);
f = (radii(1)-radii(2))/radii(2);
[lon, lat, alt] = cspice_recpgr('URANUS', ssposn, re, f);
%alt,[lon lat]*cspice_dpr

if 0
%
% compute sub-solar point
%
method = 'Near point: ellipsoid';
target = 'SUN';
fixref = 'IAU_URANUS';
%abcorr = 'LT+S';
obsrvr = 'URANUS';
[spoint, trgepc, srfvec] = cspice_subpnt(method,target,et,fixref,abcorr,obsrvr);
method = 'Near point: ellipsoid';
target = 'URANUS';
fixref = 'IAU_URANUS';
abcorr = 'LT+S';
obsrvr = 'Earth';
[ssolar, trgepc, srfvec] = cspice_subslr(method,target,et,fixeref,abcorr,obsrvr);
end

% Uranian Central Meridian Longitude
% CML is defined by the longitude of Uranus facing the Earth at a certain time.
rotate = cspice_pxform('J2000', 'IAU_URANUS', et);
sysIIIstate = rotate*state(1:3);
% modulo to get longitude
CML  = mod(atan2(sysIIIstate(2), sysIIIstate(1))*180/pi, 360);
%lat = 90 - acos(sysIIIstate(3)/sysIIIdist)*180/pi

% and Earth
target   = 'EARTH';
frame    = 'IAU_URANUS';
abcorr   = 'NONE';
observer = 'URANUS';
% in IAU_URANUS frame the rotation axis of Uranus is UranusRotAxis=(0,0,1)
[state , ltime] = cspice_spkezr(target, et, frame, abcorr, observer);

% Calculation of the angle between celestial north and Uranus rotation axis 
% as seen along the line of sight Earth-Uranus
% In IAU_URANUS coordinates system
lineOfSight = -cspice_vhat(state(1:3));
UranusRotAxis = [0;0;1];

% matrix to transform from Earth referential to Uranus referential
Earth2UranusRef = cspice_pxform('IAU_EARTH','IAU_URANUS',et);
EarthRotAxis = cspice_vhat(Earth2UranusRef*[0;0;1]);

% psi is the angle between the projections of the rotation axes of Earth and
% Uranus onto the plane perpendicular to the line of sight from Earth to
% Uranus
projUranusRotAxis = cspice_vhat( ...
                    UranusRotAxis-dot(UranusRotAxis,lineOfSight)*lineOfSight);
projEarthRotAxis = cspice_vhat( ...
                    EarthRotAxis-dot(EarthRotAxis,lineOfSight)*lineOfSight);
crossProductProj = cspice_vhat(cross(projEarthRotAxis,projUranusRotAxis));

cosa = dot(projEarthRotAxis,projUranusRotAxis);
sina = cspice_vnorm(cross(projEarthRotAxis,projUranusRotAxis));

if dot(crossProductProj,lineOfSight)<0,% projEarth,projUranus anti-clockwise
  psi = acosd(cosa);
else, % clockwise
  psi = 360-acosd(cosa);
end

seposn = state(1:3);
sedist  = norm(seposn);
% modulo to get longitude
%selong  = mod(atan2(seposn(2), seposn(1))*180/pi, 360);
selong  = atan2(seposn(2), seposn(1))*180/pi;
selat   = 90 - acos(seposn(3)/sedist)*180/pi;

% Uranian Central Meridian Longitude
% CML is defined by the longitude of Uranus facing the Earth at a certain time.
target   = 'EARTH';
frame    = 'J2000';
abcorr   = 'NONE';
observer = 'URANUS';
[state , ltime] = cspice_spkezr(target, et, frame, abcorr, observer);
rotate = cspice_pxform('J2000', 'IAU_URANUS', et);
sysIIIstate = rotate*state(1:3);
sysIIIdist = norm(sysIIIstate);
% modulo to get longitude
CML  = mod(atan2(sysIIIstate(2), sysIIIstate(1))*180/pi, 360)
%lat = 90 - acos(sysIIIstate(3)/sysIIIdist)*180/pi

% Sun-Earth in AU
sedistAU = cspice_convrt(sedist,'KM','AU');
% AU in km
AU2km = cspice_convrt(1,'AU','KM');

% Transform position of Uranus axis to Radial Tangential Normal coordinates
% Set up RTN definitions in Uranus coordinates
% R = Sun to Uranus unit vector
% T = (Omega x R) / | (Omega x R) | where Omega is Sun's spin axis  
% N completes the right-handed triad 

% normalised radial vector from Sun toward Uranus 
rvec = -cspice_vhat(ssposn);

% get the matrix that transforms position vectors from Sun to Uranus
% a specified epoch 'et'. 
% For a n-vector 'et' Sun2Uranus is an array of dimensions (3,3,n).
Sun2Uranus = cspice_pxform('IAU_SUN', 'IAU_URANUS', et);

% Sun axis orientation in Sun-centred system
sunaxis = [0.0; 0.0; 1.0];
% and in frame of Uranus
sunaxis = Sun2Uranus * sunaxis;

% Tangential is perpendicular to radial and sun axis 
tvec = cspice_vhat(cross(sunaxis, rvec));
% Normal is perpendicular to radial and tangential
nvec = cspice_vhat(cross(rvec, tvec));

% Uranus axis orientation in Uranus-centred system
jupaxis = [0.0; 0.0; 1.0]; 

% in RTN
jupaxisRTN = [rvec,tvec,nvec]'*jupaxis;
if 0
jupaxis_r = sum(jupaxis.*rvec);
jupaxis_t = sum(jupaxis.*tvec);
jupaxis_n = sum(jupaxis.*nvec);
end

% get Uranus axis / sun direction angular separation
axis_sun_ang = acos(-1.0*rvec(3))*180/pi;

if 0,
fprintf(1,'%s%s\n',char(' '*ones(1,17)), ...
        'UTC   Subsol Lat  Subsol Long   Subear Lat  Subear Long');
fprintf(1,'%s %12.6f %12.6f %12.6f %12.6f\n', ...
        cspice_et2utc(et,'C',0), sslat, sslong, selat, selong);
else
fprintf(1,'Epoch                         = %s\n', cspice_et2utc(et,'C',0));
fprintf(1,'Sub-Earth latitude, longitude = %12.4f, %12.4f\n', selat, selong);
fprintf(1,'Sub-Solar latitude, longitude = %12.4f, %12.4f\n', sslat, sslong);
end

%  It's always good form to unload kernels after use,
%  particularly in MATLAB due to data persistence.
cspice_kclear

