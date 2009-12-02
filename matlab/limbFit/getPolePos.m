function [xp,yp] = getPolePos()
% function [xp,yp] = getPolePos()

%
% $Id: getPolePos.m,v 1.2 2009/12/02 22:25:17 patrick Exp $
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

fprintf(1,'\nSelect approximate pole position and press mouse button\n');
[xp, yp] = ginput(1);

fprintf(1,'xp = %.1f yp = %.1f\n', xp, yp);
