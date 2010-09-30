function fit = selectSeeds(fit,Sx,Sy,Sls)
% function fit = selectSeeds(fit,Sx,Sy,Sls)

%
% $Id: selectSeeds.m,v 1.6 2010/09/30 18:18:13 patrick Exp $
%
% Copyright (c) 2009 
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


p     = fit.p0;

Sx    = fit.Sx;
Sy    = fit.Sy;
LS    = fit.Sls;

LSmax = fit.LSmax;
Rmin  = fit.Rmin;
Rmax  = fit.Rmax;

% find seeds with specification scale length and position
xc   = p(1);
yc   = p(2);
a    = p(3);
if length(p) == 3,
  b  = a;
  t0 = 0;
elseif length(p) == 5,
  b  = p(4);
  t0 = p(5);
end

% unit vector along major axis
xa = cosd(t0);
ya = sind(t0);

% rotate seeds along specified major axis
X =  xa*(Sx-xc) + ya*(Sy-yc);
Y = -ya*(Sx-xc) + xa*(Sy-yc);

% canonical form of ellipse (X/a)^2+(Y/b)^2=1
ellipseEq = X.^2/a^2+Y.^2/b^2;

% selection of seeds criteria 
%iSelect = find(Sls <= LSmax & ellipseEq > Rmin^2 & ellipseEq < Rmax^2);
iSelect = (Sls <= LSmax & ellipseEq > Rmin^2 & ellipseEq < Rmax^2);
fprintf(1,'Total number of seeds: %d\n', length(Sls));
fprintf(1,'Number Seeds on limb : %d\n', length(Sls(find(iSelect))))

% finally angle selection
T = 180/pi*atan2(Y, X);

if ~isempty(fit.Tlim),
  for i=1:length(fit.Tlim),
	  iSelect = iSelect & (T < fit.Tlim{i}(1) | T > fit.Tlim{i}(2) );
    fprintf(1,'Number Seeds on limb : %d\n', length(Sls(find(iSelect))));
	end
end

if 0,
% if pole position is not empty calculate the polar angle to the pole
% and substract that value to the polar angle to the seeds
if ~isempty(fit.polePos),
  Tpole = 180/pi*atan2(fit.polePos(2)./b, fit.polePos(1)./a);
	fprintf(1,'Pole polar angle %.0f deg\n', Tpole);
	T = T - Tpole;
end
if isa(fit.selectAngles,'char') | isa(fit.selectAngles,'function_handle'),
      [selectAngles, msg] = fcnchk(fit.selectAngles);
end
withinAngularSpec = selectAngles(T(iSelect));
if any(withinAngularSpec == false),
  iSelect(withinAngularSpec == false) = [];
  fprintf(1,'Number Seeds on limb : %d\n', length(iSelect))
end
end

% embed selected seeds in fit structure
fit.iSelect = find(iSelect);

