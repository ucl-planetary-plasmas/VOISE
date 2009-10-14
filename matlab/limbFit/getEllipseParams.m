function p = getEllipseParams()
% function p = getEllipseParams()

%
% $Id: getEllipseParams.m,v 1.1 2009/10/14 15:16:55 patrick Exp $
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

fprintf(1,'\nSelect approximate disc centre and press mouse button\n ');
[xc, yc] = ginput(1);

fprintf(1,'xc = %.1f yc = %.1f\n', xc, yc);

fprintf(1,'\nSelect approximate semi-major axis and press mouse button\n ');
[xa, ya] = ginput(1);

a = sqrt((xc-xa).^2+(yc-ya).^2);
ta = 180+180/pi*atan2(yc-ya,xc-xa);

fprintf(1,'a = %.1f t0 = %.0f\n', a, ta);

fprintf(1,'\nSelect approximate semi-minor axis and press mouse button\n ');
[xb, yb] = ginput(1);

b = sqrt((xc-xb).^2+(yc-yb).^2);
tb = 90+180/pi*atan2(yc-yb,xc-xb);

fprintf(1,'b = %.1f t0 = %.0f\n', b, tb);

t0 = 0.5*(ta+tb);

p = [xc,yc,a,b,t0];

