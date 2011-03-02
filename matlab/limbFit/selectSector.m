function Tlim = selectSector()
% function Tlim = selectSector()

%
% $Id: selectSector.m,v 1.1 2011/03/02 17:45:21 patrick Exp $
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

fprintf(1,'\nSelect approximate origo position and press mouse button\n');
[xo, yo] = ginput(1);

fprintf(1,'xp = %.1f yp = %.1f\n', xo, yo);

fprintf(1,'\nSelect approximate lower position sector and press mouse button\n');
[xm, ym] = ginput(1);

fprintf(1,'\nSelect approximate upper position sector and press mouse button\n');
[xM, yM] = ginput(1);

Tlim = fix(180/pi*sort([atan2(ym-yo, xm-xo), atan2(yM-yo, xM-xo)]));

fprintf(1,'Tlim = %.0f, %.0f deg\n', Tlim);
