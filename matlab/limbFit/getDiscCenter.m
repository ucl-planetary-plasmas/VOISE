function [xc,yc] = getDiscCenter()
% function [xc,yc] = getDiscCenter()

%
% $Id: getDiscCenter.m,v 1.1 2009/10/14 15:16:55 patrick Exp $
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
