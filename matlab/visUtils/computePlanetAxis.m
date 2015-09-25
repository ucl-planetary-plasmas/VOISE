function [ss,se]=computePlanetAxis(planet,epoch)
% function [ss,se]=computePlanetAxis(planet,epoch)
%
% $Id: computePlanetAxis.m,v 1.11 2015/09/25 13:45:53 patrick Exp $
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
% converts rectangular coordinates to latitudinal coordinates.
[ssrad,sslon,sslat] = cspice_reclat(ssposn);
fprintf(1,'ss   %9.f lat %+9.4f lon %+9.4f\n',ssrad,[sslat,c2grlon(sslon)]*r2d);
if 0,
ssrd1 = norm(ssposn);
ssln1 = atan2(ssposn(2), ssposn(1))*r2d;
% Get longitude in [0,360] instead of [-180,180]
ssln1 = wrap360(ssln1);
sslt1 = 90 - acos(ssposn(3)/ssrd1)*r2d;
fprintf(1,'ss1  %9.f lat %+9.4f lon %+9.4f\n',ssrd1,[sslt1,c2gdlon(ssln1)]);
fprintf(1,'diff %+9.2g     %+9.2g     %+9.2g\n',...
        ssrad-ssrd1,sslat*r2d-sslt1,sslon*r2d-ssln1)
end
% converts rectangular coordinates to planetographic coordinates
[lon,lat,alt] = cspice_recpgr(PLANET,ssposn,re,f);
fprintf(1,'alt  %9.f lat %+9.4f lon %+9.4f\n',alt,[lat,lon]*r2d);

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

fprintf(1,'spc  %9.f lat %+9.4f lon %+9.4f\n',spcrad,[spclat,c2grlon(spclon)]*r2d);
fprintf(1,'spg  %9.f lat %+9.4f lon %+9.4f\n',spgalt,[spglat,spglon]*r2d);
fprintf(1,'supc %9.f lat %+9.4f lon %+9.4f\n',supcrd,[supclt,c2grlon(supcln)]*r2d);
fprintf(1,'supg %9.f lat %+9.4f lon %+9.4f\n',supgal,[supglt,supgln]*r2d);

% Return sub-solar point parameters
[ss.rad,ss.lon,ss.lat] = deal(ssrad,sslon*r2d,sslat*r2d); 
[ss.rad,ss.lon,ss.lat] = deal(supcrd,c2grlon(supcln)*r2d,supclt*r2d); 
ss.trgepc = trgepc;

% and Sub-Earth Point
fprintf(1,'Sub-Earth Point\n');
% Get position of Earth with respect to planet 
target   = 'EARTH';
obsrvr   = PLANET;
[state,ltime] = cspice_spkezr(target,et,fixref,abcorr,obsrvr);
seposn = state(1:3);
[serad,selon,selat] = cspice_reclat(seposn);
fprintf(1,'se   %9.f lat %+9.4f lon %+9.4f\n',serad,[selat,selon]*r2d);
if 0,
serd1  = norm(seposn);
seln1  = atan2(seposn(2), seposn(1))*r2d;
% Get longitude in [0,360] instead of [-180,180]
seln1  = wrap360(seln1);
selt1   = 90 - acos(seposn(3)/serd1)*r2d;
fprintf(1,'se1  %9.f lat %+9.4f lon %+9.4f\n',serd1,[selt1,seln1]);
fprintf(1,'diff %+9.2g     %+9.2g     %+9.2g\n',...
        serad-serd1,selat*r2d-selt1,selon*r2d-seln1)
end

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

fprintf(1,'spc  %9.f lat %+9.4f lon %+9.4f\n',spcrad,[spclat,c2grlon(spclon)]*r2d);
fprintf(1,'spg  %9.f lat %+9.4f lon %+9.4f\n',spgalt,[spglat,spglon]*r2d);
fprintf(1,'opc  %9.f lat %+9.4f lon %+9.4f\n',opcrad,[opclat,c2grlon(opclon)]*r2d);
fprintf(1,'opg  %9.f lat %+9.4f lon %+9.4f\n',opgalt,[opglat,opglon]*r2d);

% Return sub-Earth point parameters
[se.rad,se.lon,se.lat] = deal(serad,selon*r2d,selat*r2d); 
[se.rad,se.lon,se.lat] = deal(opcrad,c2grlon(opclon)*r2d,opclat*r2d); 
se.trgepc = trgepc;

% Planet-Earth in AU
se.distkm = opcrad;
se.distAU = cspice_convrt(se.distkm,'KM','AU');
% AU in km
AU2km = cspice_convrt(1,'AU','KM');

% Planet's Central Meridian Longitude
% Longitude of the planet facing the Earth at a certain time
if 0,
% Get position of Earth with respect to planet
target   = 'EARTH';
frame    = 'J2000';
obsrvr   = PLANET;
[state,ltime] = cspice_spkezr(target,et,frame,abcorr,obsrvr);
rotate = cspice_pxfrm2('J2000',IAU_PLANET,et-ltime,et);
sysIIIstate = rotate*state(1:3);
[~,CML1,~] = cspice_reclat(sysIIIstate(1:3));
fprintf(1,'CML(III)                 = %+9.4f deg\n', CML1*r2d);
end
CML = c2grlon(opclon);
fprintf(1,'CML(III)                 = %+9.4f deg\n', CML*r2d);

% Return CML
se.CML = CML;


% Calculation of the apparent angle between celestial North and the planet's
% rotation axis seen from Earth 

% nSky: normal vector to Earth sky plane pointing from planet to Earth
if 0,
% seposn: position of Earth with respect to planet
nSky = cspice_vhat(seposn);
% matrix to transform from Earth referential to planet's referential
Earth2PlanetRef = cspice_pxfrm2('IAU_EARTH',IAU_PLANET,et,et);
else
% obspos: observer's position relative to the center of the target
nSky = cspice_vhat(obspos);
% matrix to transform from Earth referential to planet's referential
Earth2PlanetRef = cspice_pxfrm2('IAU_EARTH',IAU_PLANET,et,se.trgepc);
end

% Planet's rotation axis 
planetRotAxis = [0;0;1];

% normalised direction of Earth rotation axis 
EarthRotAxis = cspice_vhat(Earth2PlanetRef*[0;0;1]);

% psi is the oriented angle (positive anti-clockwise) between the projections
% of Earth and the planet's rotation axes onto the Earth sky plane 
projPlanetRotAxis = projVecOnPlane(planetRotAxis,nSky);
projEarthRotAxis = projVecOnPlane(EarthRotAxis,nSky);

% cosine and sine of psi in a frame with anti-clockwise orientation
cosa = dot(projEarthRotAxis,projPlanetRotAxis);
sina = dot(cross(projEarthRotAxis,projPlanetRotAxis),nSky);
psi = atan2(sina,cosa);
fprintf(1,'psi(NP)                  = %+9.4f deg\n', psi*r2d);

if 0,
if sina>0, % Angle(projEarthRotAxis,projPlanetRotAxis)>0 psi in [0,pi] 
  psi1 = acos(cosa);
else, % Angle(projEarthRotAxis,projPlanetRotAxis)<0 psi in [-pi,0]
  psi1 = -acos(cosa);
end
fprintf(1,'psi1(NP)                 = %+9.4f deg\n', psi1*r2d);
fprintf(1,'diff                     = %+9.2f deg\n', (psi-psi1)*r2d);
end

% Return psi
se.psi = psi;

% Transform position of planet's axis to Radial Tangential Normal coordinates
% Set up RTN definitions in planet's coordinates
% R = Sun to planet unit vector
% T = (Omega x R) / | (Omega x R) | where Omega is Sun's spin axis  
% N completes the right-handed triad 

% normalised radial vector from Sun toward planet 
rvec = -cspice_vhat(ssposn);

% get the matrix that transforms position from Sun at epoch from to planet at epoch to 
% For n-vectors from/to Sun2Planet is an array of dimensions (3,3,n).
Sun2Planet = cspice_pxfrm2('IAU_SUN',IAU_PLANET,ss.trgepc,ss.trgepc);

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
if 0,
planetAxis_r = sum(planetAxis.*rvec);
planetAxis_t = sum(planetAxis.*tvec);
planetAxis_n = sum(planetAxis.*nvec);
end

% get planet's axis / sun direction angular separation
planetAxisSunAng = acos(-1.0*rvec(3));

fprintf(1,'planetAxisSunAng         = %+9.4f deg\n',planetAxisSunAng*r2d);
fprintf(1,'90-ss.lat                = %+9.4f deg\n',90-ss.lat);

fprintf(1,'Earth epoch              = %s\n',cspice_et2utc(et,'C',0));
fprintf(1,'Target epoch             = %s\n',cspice_et2utc(se.trgepc,'C',0));
fprintf(1,'Sub-Earth lat, lon (deg) = %+9.4f, %+9.4f\n',se.lat,se.lon);
fprintf(1,'Sub-Solar lat, lon (deg) = %+9.4f, %+9.4f\n',ss.lat,ss.lon);
fprintf(1,'Earth-Sun lon difference = %+9.4f deg\n',ss.lon-se.lon);
fprintf(1,'Earth-%s distance%s= %9.4f AU (%.0f km)\n',...
        Planet,char(' '*ones(1,10-length(Planet))),se.distAU,se.distkm);

%  It's always good form to unload kernels after use,
%  particularly in MATLAB due to data persistence.
cspice_kclear

% Get longitude in [0,360] instead of [-180,180]
function lon = wrap360(lon)
lon  = mod(lon, 360);

% Convert planetocentric to planetographic longitude in rad
function lon = c2grlon(lon)
lon = mod(2*pi - lon, 2*pi);

% Convert planetocentric to planetographic longitude in deg
function lon = c2gdlon(lon)
lon = mod(360 - lon, 360);

% Compute unit vector of 
% the projection of x onto the plane defined by its normal vector np
function xproj = projVecOnPlane(x,np)
% substract to x the quantity dot(x,np)*np, the component of x along np
xproj = cspice_vhat(x - dot(x,np)*np);

