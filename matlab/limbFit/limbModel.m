function r = limbModel(t,p)
% function r = limbModel(t,p)

%
% $Id: limbModel.m,v 1.1 2009/10/16 13:55:35 patrick Exp $
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

if length(p) == 3,

  r = circle(t,p);

elseif length(p) == 5,

  r = ellipse(t,p);

end
