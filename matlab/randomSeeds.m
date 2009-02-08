function S = randomSeeds(nr,nc,ns,pc)
% function S = randomSeeds(nr,nc,ns,pc)

%
% $Id: randomSeeds.m,v 1.1 2009/02/08 21:07:16 patrick Exp $
%
% Copyright (c) 2008 
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

% initialise array S(ns,2) 
% seed s has coordinates (x,y) = S(s, 1:2) 

if ~exist('pc','var'),
	pc = 0.02;
end

% no seeds in the 100*pc % from the boundary of the image
S = round([(nc-2*pc*nc)*(rand(ns,1))+pc*nc, ...
           (nr-2*pc*nr)*(rand(ns,1))+pc*nr]);


