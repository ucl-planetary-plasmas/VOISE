function [a,b,e] = getPlanetGeometry(planetName)

% generic kernel path
persistent spiceKernelsPath

spiceKernelsPath = '/home/patrick/research/codes/spice/data/';

% Load orientation data for planets, natural 
% satellites, the Sun, and selected asteroids
% ftp://naif.jpl.nasa.gov/pub/naif/generic_kernels/pck/
cspice_furnsh([spiceKernelsPath 'pck00008.tpc']);

radii = cspice_bodvrd(planetName,'RADII',3);

a = radii(1);
b = radii(3);
e = sqrt(a^2-b^2)/a;

%  It's always good form to unload kernels after use,
%  particularly in MATLAB due to data persistence.
cspice_kclear


