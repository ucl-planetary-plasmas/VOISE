function [ss,se,CML,psi]=computePlanetAxis(planet,epoch)
% function [ss,se,CML,psi]=computePlanetAxis(planet,epoch)
%
% $Id: computePlanetAxis.m,v 1.8 2015/09/23 17:49:29 patrick Exp $
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
r2d = cspice_dpr;

% Correct for one-way light time and stellar aberration using a Newtonian
% formulation. This option modifies the state obtained with the 'LT' option to
% account for the observer's velocity relative to the solar system barycenter.
% The result is the apparent sub-observer point as seen by the observer
%abcorr   = 'LT+S';

% Converged Newtonian light time and stellar aberration corrections. This
% option produces a solution that is at least as accurate at that obtainable
% with the 'LT+S' option. Whether the 'CN+S' solution is substantially more
% accurate depends on the geometry of the participating objects and on the
% accuracy of the input data. In all cases this routine will execute more
% slowly when a converged solution is computed
abcorr   = 'CN+S';

% body-fixed, body-centered reference frame associated with the target body
fixref = IAU_PLANET;

% Nearest point on the target relative to the Sun/observer
%method = 'Near point: ellipsoid';

% Target surface intercept of the line containing the Sun/observer and
% the target's center
method = 'Intercept:  ellipsoid';

% load necessary kernels
loadPlanetSpiceKernels(planet);

% re, f are equatorial radius and flattening of the planet
%  necessary for cspice_recpgr()
radii = cspice_bodvrd(PLANET,'RADII',3);
re = radii(1);
f = (radii(1)-radii(3))/radii(1);

% convert UTC to ephemeris time 
et = cspice_str2et(epoch);

% Sub-solar point calculations
fprintf(1,'Sub-Solar Point\n');

% Get position of the Sun with respect to planet 
target = 'SUN'; 
obsrvr = PLANET;
% state=[x,y,z,vx,vy,vz]' (6xn) and light time 'ltime' n-vector 
% for ephemeris time 'et' n-vector.
target = 'SUN';
obsrvr = PLANET;
[state,ltime] = cspice_spkezr(target,et,fixref,abcorr,obsrvr); 
% converts rectangular coordinates to latitudinal coordinates
ssposn = state(1:3);
ssdist1 = norm(ssposn);
sslon1 = atan2(ssposn(2), ssposn(1))*r2d;
% Get longitude in [0,360] instead of [-180,180]
%sslon1  = wrap360(sslon1);
sslat1 = 90 - acos(ssposn(3)/ssdist1)*r2d;
fprintf(1,'ss1  %10.f lat %+9.5f lon %+9.5f\n',ssdist1,[sslat1,sslon1]);
% converts rectangular coordinates to latitudinal coordinates.
[ssdist,sslon,sslat] = cspice_reclat(ssposn);
fprintf(1,'ss   %10.f lat %+9.5f lon %+9.5f\n',ssdist,[sslat,c2grlon(sslon)]*r2d);
%fprintf(1,'delta  %+10g  %+12.5g   %+12.5g\n',...
%        ssdist-ssdist1,sslat*r2d-sslat1,sslon*r2d-sslon1)
% converts rectangular coordinates to planetographic coordinates
[lon,lat,alt] = cspice_recpgr(PLANET,ssposn,re,f);
fprintf(1,'alt  %10.f lat %+9.5f lon %+9.5f\n',alt,[lat,lon]*r2d);

% Compute sub-solar point with spice function
target = PLANET;
obsrvr = 'EARTH';
[spoint,trgepc,srfvec] = cspice_subslr(method,target,et,fixref,abcorr,obsrvr);
% Convert sub-solar point's rectangular coordinates to
% planetocentric radius, longitude and latitude
[spcrad,spclon,spclat] = cspice_reclat(spoint);
% and to planetographic longitude, latitude and altitude
[spglon,spglat,spgalt] = cspice_recpgr(target,spoint,re,f);
% Compute Sun's apparent position relative to the target's center at trgepc
[sunpos,sunlt] = cspice_spkpos('SUN',trgepc,fixref,abcorr,target);
% Express the Sun's location in planetocentric coordinates
[supcrd,supcln,supclt] = cspice_reclat(sunpos);
% and to planetographic longitude, latitude and altitude
[supgln,supglt,supgal] = cspice_recpgr(target,sunpos,re,f);

fprintf(1,'spc  %10.f lat %+9.5f lon %+9.5f\n',spcrad,[spclat,c2grlon(spclon)]*r2d);
fprintf(1,'spg  %10.f lat %+9.5f lon %+9.5f\n',spgalt,[spglat,spglon]*r2d);
fprintf(1,'supc %10.f lat %+9.5f lon %+9.5f\n',supcrd,[supclt,c2grlon(supcln)]*r2d);
fprintf(1,'supg %10.f lat %+9.5f lon %+9.5f\n',supgal,[supglt,supgln]*r2d);

% Return sub-solar point parameters
[ss.lon,ss.lat] = deal(sslon*r2d,sslat*r2d); 
%[ss.lon,ss.lat] = deal(c2grlon(supcln)*r2d,supclt*r2d); 

% and Sub-Earth Point
fprintf(1,'Sub-Earth Point\n');
% Get position of Earth with respect to planet 
target   = 'EARTH';
obsrvr   = PLANET;
[state,ltime] = cspice_spkezr(target,et,fixref,abcorr,obsrvr);
seposn = state(1:3);
sedist1  = norm(seposn);
selon1  = atan2(seposn(2), seposn(1))*r2d;
% Get longitude in [0,360] instead of [-180,180]
%selon1  = wrap360(selon1);
selat1   = 90 - acos(seposn(3)/sedist1)*r2d;
%fprintf(1,'se1  %10.f lat %+9.5f lon %+9.5f\n',sedist1,[selat1,selon1]);
[sedist,selon,selat] = cspice_reclat(seposn);
fprintf(1,'se   %10.f lat %+9.5f lon %+9.5f\n',sedist,[selat,selon]*r2d);
%fprintf(1,'delta  %+10g  %+12.5g   %+12.5g\n',...
%        sedist-sedist1,selat*r2d-selat1,selon*r2d-selon1)

target = PLANET;
obsrvr = 'EARTH';
% Compute the sub-observer point
[spoint,trgepc,srfvec] = cspice_subpnt(method,target,et,fixref,abcorr,obsrvr);

% Compute the observer's distance from SPOINT
odist = norm(srfvec);
% Convert sub-observer point's rectangular coordinates to
% planetocentric radius, longitude and latitude
[spcrad,spclon,spclat] = cspice_reclat(spoint);
% and to planetographic longitude, latitude and altitude
[spglon,spglat,spgalt] = cspice_recpgr(target,spoint,re,f);
% Compute the observer's position relative to the center of the
% target, where the center's location has been adjusted using
% the aberration corrections applicable to the sub-point.
obspos = spoint - srfvec;
% Express the observer's location in planetocentric coordinates
[opcrad,opclon,opclat] = cspice_reclat(obspos);
% and to planetographic longitude, latitude and altitude
[opglon,opglat,opgalt] = cspice_recpgr(target,obspos,re,f);

fprintf(1,'spc  %10.f lat %+9.5f lon %+9.5f\n',spcrad,[spclat,c2grlon(spclon)]*r2d);
fprintf(1,'spg  %10.f lat %+9.5f lon %+9.5f\n',spgalt,[spglat,spglon]*r2d);
fprintf(1,'opc  %10.f lat %+9.5f lon %+9.5f\n',opcrad,[opclat,c2grlon(opclon)]*r2d);
fprintf(1,'opg  %10.f lat %+9.5f lon %+9.5f\n',opgalt,[opglat,opglon]*r2d);

% Return sub-Earth point parameters
[se.lon,se.lat] = deal(selon*r2d,selat*r2d); 
[se.lon,se.lat] = deal(c2grlon(opclon)*r2d,opclat*r2d); 

% Planet-Earth in AU
se.distkm = opcrad;
se.distAU = cspice_convrt(sedist,'KM','AU');
% AU in km
AU2km = cspice_convrt(1,'AU','KM');

% planet's Central Meridian Longitude
% CML defined as longitude of the planet facing the Earth at a certain time
target   = 'EARTH';
frame    = 'J2000';
obsrvr   = PLANET;
[state,ltime] = cspice_spkezr(target,et,frame,abcorr,obsrvr);
rotate = cspice_pxform('J2000',IAU_PLANET,et);
sysIIIstate = rotate*state(1:3);
CML1 = atan2(sysIIIstate(2), sysIIIstate(1))*r2d;
% Get longitude in [0,360] instead of [-180,180]
CML1 = wrap360(CML1); 
[~,CML,~] = cspice_reclat(sysIIIstate(1:3));
%[CML,~,~] = cspice_recpgr(obsrvr,sysIIIstate(1:3),re,f);
CML = c2grlon(opclon);
fprintf(1,'CML(III) %+9.5f %+9.5f deg\n', CML1, CML*r2d);


% Calculation in planet's IAU coordinates of the angle between celestial North
% and planet's rotation axis as seen along the line of sight Earth to planet
% Normalised pointing direction
lineOfSight = -cspice_vhat(seposn);
% Planet's rotation axis 
planetRotAxis = [0;0;1];

% matrix to transform from Earth referential to planet's referential
Earth2PlanetRef = cspice_pxform('IAU_EARTH',IAU_PLANET,et);
% Normalised direction of Earth rotation axis 
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
  psi = acos(cosa);
else, % clockwise
  psi = 360-acos(cosa);
end
fprintf(1,'psi              %+9.5f deg\n', psi*r2d);


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
sunaxis = [0;0;1];
% and in frame of planet
sunaxis = Sun2Planet * sunaxis;

% Tangential is perpendicular to radial and sun axis 
tvec = cspice_vhat(cross(sunaxis, rvec));
% Normal is perpendicular to radial and tangential
nvec = cspice_vhat(cross(rvec, tvec));

% Planet's axis orientation in planet-centred system
planetAxis = [0;0;1]; 

% in RTN
planetAxisRTN = [rvec,tvec,nvec]'*planetAxis;
if 0
planetAxis_r = sum(planetAxis.*rvec);
planetAxis_t = sum(planetAxis.*tvec);
planetAxis_n = sum(planetAxis.*nvec);
end

% get planet's axis / sun direction angular separation
planetAxisSunAng = acos(-1.0*rvec(3));

fprintf(1,'planetAxisSunAng %+9.5f deg\n', planetAxisSunAng*r2d);

if 0,
fprintf(1,'%s%s\n',char(' '*ones(1,17)), ...
        'UTC   Subsol Lat  Subsol Long   Subear Lat  Subear Long');
fprintf(1,'%s %12.6f %12.6f %12.6f %12.6f\n', ...
        cspice_et2utc(et,'C',0), ss.lat, ss.lon, se.lat, se.lon);
else
fprintf(1,'Epoch                         = %s\n', cspice_et2utc(et,'C',0));
fprintf(1,'Sub-Earth latitude, longitude = %12.4f, %12.4f\n', se.lat, se.lon);
fprintf(1,'Sub-Solar latitude, longitude = %12.4f, %12.4f\n', ss.lat, ss.lon);
fprintf(1,'Earth-%s distance = %.4f AU (%.0f km)\n',Planet,se.distAU,se.distkm);
end

%  It's always good form to unload kernels after use,
%  particularly in MATLAB due to data persistence.
cspice_kclear

% Get longitude in [0,360] instead of [-180,180]
function lon = wrap360(lon)
lon  = mod(lon, 360);

% Convert planetocentric to planetographic longitude
function lon = c2grlon(lon)
lon = mod(2*pi - lon, 2*pi);
