function flag = ex2Angles(T)
% function flag = ex2Angles(T)

%
% $Id: ex2Angles.m,v 1.1 2009/11/06 17:15:23 patrick Exp $
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

% do not consider points at the (South) pole
flag = T<=-90-20 | T>=-90+20;

fprintf('exAngles: Removed %d seeds at the pole\n', length(find(flag==false)));
