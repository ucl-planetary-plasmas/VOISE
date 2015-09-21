function [ss,se,CML,psi,sedistAU,AU2km]=computePlanetAxis(planet,epoch)
% function [ss,se,CML,psi,sedistAU,AU2km]=computePlanetAxis(planet,epoch)
%
% $Id: computePlanetAxis.m,v 1.7 2015/09/21 14:03:55 patrick Exp $
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

planet = lower(planet);
PLANET = upper(planet);
Planet = [PLANET(1) planet(2:end)];
IAU_PLANET = ['IAU_' PLANET];

% conversion factor
deg2rad = cspice_dpr;

%abcorr   = 'LT';
%abcorr   = 'LT+S';
abcorr   = 'CN+S';

method = 'Near point: ellipsoid';
%method = 'Intercept:  ellipsoid';

% load necessary kernels
loadPlanetSpiceKernels(planet);

% convert UTC to ephemeris time 
et = cspice_str2et(epoch);

% Sub-solar point calculation
% Get position of Sun with respect to planet 
fprintf(1,'Sub-Solar Point\n');
target = 'SUN';
frame  = IAU_PLANET;
obsrvr = PLANET;
% state=[x,y,z,vx,vy,vz]' matrix (6xn) and light time 'ltime' n-vector 
% for ephemeris time 'et' n-vector.
[state,ltime] = cspice_spkezr(target,et,frame,abcorr,obsrvr);
%[state,ltime] = cspice_spkezr(obsrvr,et,frame,abcorr,target);
% converts rectangular coordinates to latitudinal coordinates
ssposn = state(1:3);
%ssdist1 = norm(ssposn);
%ss.lon1 = atan2(ssposn(2), ssposn(1))*deg2rad;
%ss.lat1 = 90 - acos(ssposn(3)/ssdist1)*deg2rad;
%fprintf(1,'ssdist %10.f lat %+9.5f lon %+9.5f\n',ssdist1,[ss.lat1,ss.lon1]);
% converts rectangular coordinates to latitudinal coordinates.
[ssdist,ss.lon,ss.lat] = cspice_reclat(ssposn);
[ss.lon,ss.lat] = deal(ss.lon*deg2rad,ss.lat*deg2rad);
fprintf(1,'ss   %10.f lat %+9.5f lon %+9.5f\n',ssdist,[ss.lat,ss.lon]);
%fprintf(1,'delta  %+10g  %+12.5g   %+12.5g\n',ssdist-ssdist1,ss.lat-ss.lat1,ss.lon-ss.lon1)
% re, f are equatorial radius and flattening of the planet
radii = cspice_bodvrd(PLANET,'RADII',3);
re = radii(1);
f = (radii(1)-radii(3))/radii(1);
% converts rectangular coordinates to planetographic coordinates
[lon,lat,alt] = cspice_recpgr(PLANET,ssposn,re,f);
%fprintf(1,'alt    %10.f lat %+9.5f lon %+9.5f\n',alt,[lat,lon]*deg2rad);

target = PLANET;
fixref = IAU_PLANET;
obsrvr = 'EARTH';
% Compute sub-solar point using light time and stellar aberration corrections
[spoint,trgepc,srfvec] = cspice_subslr(method,target,et,fixref,abcorr,obsrvr);

[spcrad,spclon,spclat] = cspice_reclat(spoint);
[spglon,spglat,spgalt] = cspice_recpgr(target,spoint,re,f);

[sunpos,sunlt] = cspice_spkpos('sun',trgepc,fixref,abcorr,target);

[supcrd,supcln,supclt] = cspice_reclat(sunpos);
[supgln,supglt,supgal] = cspice_recpgr(target,sunpos,re,f);

fprintf(1,'spc  %10.f lat %+9.5f lon %+9.5f\n',spcrad,[spclat,spclon]*deg2rad);
fprintf(1,'spg  %10.f lat %+9.5f lon %+9.5f\n',spgalt,[spglat,spglon]*deg2rad);
fprintf(1,'supc %10.f lat %+9.5f lon %+9.5f\n',supcrd,[supclt,supcln]*deg2rad);
fprintf(1,'supg %10.f lat %+9.5f lon %+9.5f\n',supgal,[supglt,supgln]*deg2rad);

% and Sub Earth Point
fprintf(1,'Sub-Earth Point\n');
% Get position of Earth with respect to planet 
target   = 'EARTH';
frame    = IAU_PLANET;
obsrvr   = PLANET;
% in planet's IAU frame the rotation axis of the planet is planetRotAxis=(0,0,1)
[state,ltime] = cspice_spkezr(target,et,frame,abcorr,obsrvr);

seposn = state(1:3);
sedist  = norm(seposn);
% modulo to get longitude
%se.lon  = mod(atan2(seposn(2), seposn(1))*deg2rad, 360);
se.lon  = atan2(seposn(2), seposn(1))*deg2rad;
se.lat   = 90 - acos(seposn(3)/sedist)*deg2rad;
fprintf(1,'se   %10.f lat %+9.5f lon %+9.5f\n',sedist,[se.lat,se.lon]);

target = PLANET;
fixref = IAU_PLANET;
obsrvr = 'EARTH';
[spoint,trgepc,srfvec] = cspice_subpnt(method,target,et,fixref,abcorr,obsrvr);

[spcrad,spclon,spclat] = cspice_reclat(spoint);
[spglon,spglat,spgalt] = cspice_recpgr(target,spoint,re,f);

obspos = spoint - srfvec;
[opcrad,opclon,opclat] = cspice_reclat(obspos);
[opglon,opglat,opgalt] = cspice_recpgr(target,obspos,re,f);

fprintf(1,'spc  %10.f lat %+9.5f lon %+9.5f\n',spcrad,[spclat,spclon]*deg2rad);
fprintf(1,'spg  %10.f lat %+9.5f lon %+9.5f\n',spgalt,[spglat,spglon]*deg2rad);
fprintf(1,'opc  %10.f lat %+9.5f lon %+9.5f\n',opcrad,[opclat,opclon]*deg2rad);
fprintf(1,'opg  %10.f lat %+9.5f lon %+9.5f\n',opgalt,[opglat,opglon]*deg2rad);

% Calculation of the angle between celestial North and planet's rotation axis 
% as seen along the line of sight from Earth to planet
% In planet's IAU coordinates system
lineOfSight = -cspice_vhat(seposn);
planetRotAxis = [0;0;1];

% matrix to transform from Earth referential to planet's referential
Earth2PlanetRef = cspice_pxform('IAU_EARTH',IAU_PLANET,et);
EarthRotAxis = cspice_vhat(Earth2PlanetRef*[0;0;1]);

% psi is the angle between the projections of the rotation axes of Earth and
% planet onto the plane perpendicular to the line of sight from Earth to
% the planet
projPlanetRotAxis = cspice_vhat( ...
                    planetRotAxis-dot(planetRotAxis,lineOfSight)*lineOfSight);
projEarthRotAxis = cspice_vhat( ...
                    EarthRotAxis-dot(EarthRotAxis,lineOfSight)*lineOfSight);
crossProductProj = cspice_vhat(cross(projEarthRotAxis,projPlanetRotAxis));

cosa = dot(projEarthRotAxis,projPlanetRotAxis);
sina = cspice_vnorm(cross(projEarthRotAxis,projPlanetRotAxis));

if dot(crossProductProj,lineOfSight)<0,% projEarth,projPlanet anti-clockwise
  psi = acosd(cosa);
else, % clockwise
  psi = 360-acosd(cosa);
end
fprintf(1,'psi       %+9.5f deg\n', psi);


%%%%% WORKING!!!!
% planet's Central Meridian Longitude
% CML is defined as longitude of the planet facing the Earth at a certain time
target   = 'EARTH';
frame    = 'J2000';
obsrvr   = PLANET;
[state,ltime] = cspice_spkezr(target,et,frame,abcorr,obsrvr);
rotate = cspice_pxform('J2000',IAU_PLANET,et);
sysIIIstate = rotate*state(1:3);
sysIIIdist = norm(sysIIIstate);
% modulo to get longitude
CML  = mod(atan2(sysIIIstate(2), sysIIIstate(1))*deg2rad, 360);
fprintf(1,'CML(III) %+9.5f deg\n', CML);
%lat = 90 - acos(sysIIIstate(3)/sysIIIdist)*deg2rad
%%%%%% WORKING!!!!

% Planet-Earth in AU
sedistAU = cspice_convrt(sedist,'KM','AU');
% AU in km
AU2km = cspice_convrt(1,'AU','KM');

% Transform position of planet's axis to Radial Tangential Normal coordinates
% Set up RTN definitions in planet's coordinates
% R = Sun to planet unit vector
% T = (Omega x R) / | (Omega x R) | where Omega is Sun's spin axis  
% N completes the right-handed triad 

% normalised radial vector from Sun toward planet 
rvec = -cspice_vhat(ssposn);

% get the matrix that transforms position vectors from Sun to planet
% a specified epoch 'et'. 
% For a n-vector 'et' Sun2Planet is an array of dimensions (3,3,n).
Sun2Planet = cspice_pxform('IAU_SUN', IAU_PLANET, et);

% Sun axis orientation in Sun-centred system
sunaxis = [0.0; 0.0; 1.0];
% and in frame of planet
sunaxis = Sun2Planet * sunaxis;

% Tangential is perpendicular to radial and sun axis 
tvec = cspice_vhat(cross(sunaxis, rvec));
% Normal is perpendicular to radial and tangential
nvec = cspice_vhat(cross(rvec, tvec));

% Planet's axis orientation in planet-centred system
jupaxis = [0.0; 0.0; 1.0]; 

% in RTN
jupaxisRTN = [rvec,tvec,nvec]'*jupaxis;
if 0
jupaxis_r = sum(jupaxis.*rvec);
jupaxis_t = sum(jupaxis.*tvec);
jupaxis_n = sum(jupaxis.*nvec);
end

% get planet's axis / sun direction angular separation
axis_sun_ang = acos(-1.0*rvec(3))*deg2rad;

if 0,
fprintf(1,'%s%s\n',char(' '*ones(1,17)), ...
        'UTC   Subsol Lat  Subsol Long   Subear Lat  Subear Long');
fprintf(1,'%s %12.6f %12.6f %12.6f %12.6f\n', ...
        cspice_et2utc(et,'C',0), ss.lat, ss.lon, se.lat, se.lon);
else
fprintf(1,'Epoch                         = %s\n', cspice_et2utc(et,'C',0));
fprintf(1,'Sub-Earth latitude, longitude = %12.4f, %12.4f\n', se.lat, se.lon);
fprintf(1,'Sub-Solar latitude, longitude = %12.4f, %12.4f\n', ss.lat, ss.lon);
fprintf(1,'Earth-%s distance = %.4f AU (%.0f km)\n',Planet,sedistAU,sedistAU*AU2km);
end

%  It's always good form to unload kernels after use,
%  particularly in MATLAB due to data persistence.
cspice_kclear

