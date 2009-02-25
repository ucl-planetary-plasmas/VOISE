function [sslat, sslong, selat, selong] = computeJupiterAxis(UDATE)

% Load a leapseconds kernel.
% ftp://naif.jpl.nasa.gov/pub/naif/generic_kernels/lsk/
cspice_furnsh('naif0009.tls');

% Load planetary ephemeris 
% ftp://naif.jpl.nasa.gov/pub/naif/generic_kernels/spk/planets/
cspice_furnsh('de421.bsp');

% Load satellite ephemeris 
% ftp://naif.jpl.nasa.gov/pub/naif/generic_kernels/spk/satellites/
cspice_furnsh('jup263.bsp');

% Load orientation data for planets, natural 
% satellites, the Sun, and selected asteroids
% ftp://naif.jpl.nasa.gov/pub/naif/generic_kernels/pck/
cspice_furnsh('pck00008.tpc');


% Convert string date to cspice
et = cspice_str2et(UDATE);


% Get position of Sun with respect to Jupiter
target   = 'SUN';
frame    = 'IAU_JUPITER';
abcorr   = 'NONE';
observer = 'JUPITER';
% Look up the 'state' vectors and light time values
% 'ltime'  corresponding to the vector of input
% ephemeris time 'et'.
[state , ltime] = cspice_spkezr(target, et, frame, abcorr, observer);

ssposn = state(1:3)';
ssdist  = norm(ssposn);
sslong  = atan2(ssposn(2), ssposn(1))*180/pi;
sslat = 90 - acos(ssposn(3)/ssdist)*180/pi;
if sslong < 0,
  sslong = sslong + 360;
end

% and Earth
target   = 'EARTH';
frame    = 'IAU_JUPITER';
abcorr   = 'NONE';
observer = 'JUPITER';
[state , ltime] = cspice_spkezr(target, et, frame, abcorr, observer);

seposn = state(1:3);
sedist  = norm(seposn);
selong  = atan2(seposn(2), seposn(1))*180/pi;
selat   = 90 - acos(seposn(3)/sedist)*180/pi;
if selong < 0,
  selong = selong + 360;
end

% Transform position of Jupiter axis to RTN coordinates
% Set up RTN definitions in Jupiter coordinates

% R vector is normalized (and -ive) version of ssposn above
rvec = -1.0* (ssposn / ssdist);

% Convert the n-vector of 'et' to an array of corresponding
% transformation matrices (dimensions (3,3,n) ).
XForm = cspice_pxform('IAU_SUN', 'IAU_JUPITER', et);

%sunaxis = cspice_mxv(XForm,[0.0, 0.0, 1.0])
sunaxis = (XForm * [0.0, 0.0, 1.0]')';

%tvec = cspice_vcrss(sunaxis, rvec);
tvec = cross(sunaxis, rvec);
tvec = cspice_vhat(tvec')';

%nvec = cspice_vcrss(rvec, tvec);
nvec = cross(rvec, tvec);
nvec = cspice_vhat(nvec')';


jupaxis = [0.0, 0.0, 1.0]; % Jupiter axis orientation in Jupiter-centred system

jupaxis_r = sum(jupaxis.*rvec);
jupaxis_t = sum(jupaxis.*tvec);
jupaxis_n = sum(jupaxis.*nvec);

% get Jupiter axis / sun dirn angular separation
axis_sun_ang = acos(-1.0*rvec(3))*180/pi;

fprintf(1,'%s%s\n',char(' '*ones(1,16)), ...
        'UTC  Subsol Lat Subsol Long Subear Lat Subear Long');
fprintf(1,'%s %12.6f %12.6f %12.6f %12.6f\n', ...
        UDATE, sslat, sslong, selat, selong);

%  It's always good form to unload kernels after use,
%  particularly in MATLAB due to data persistence.
cspice_kclear

