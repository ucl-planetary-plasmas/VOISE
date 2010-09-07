function dr = jacLimbModel(t,r,p,dp,func)
% function dr = jacLimbModel(t,r,p,dp,func)

%
% $Id: jacLimbModel.m,v 1.2 2010/09/07 17:58:42 patrick Exp $
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

  dr = dcircle(t,r,p,dp,func);

elseif length(p) == 5,

  dr = dellipse(t,r,p,dp,func);

end

