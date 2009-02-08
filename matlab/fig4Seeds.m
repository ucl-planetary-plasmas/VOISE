function S = fig4Seeds(nr,nc,ns,pc)
% function S = fig4Seeds(nr,nc,ns,pc)

%
% $Id: fig4Seeds.m,v 1.1 2009/02/08 21:07:15 patrick Exp $
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

ns = 9;

x = round(nc*[1/4; 1/2; 3/4]);
y = round(nr*[1/4; 1/2; 3/4]);

S = [[x(1:3); x([1 3]); x(1:3); x(2)], ...
     [y(3)*ones(3,1); y(2)*ones(2,1); y(1)*ones(3,1); y(2)]];

if ~exist('pc','var'),
  pc = 0.25;
end

if pc, % random fluctuation of 100*pc % of distance between seeds
  r = round([pc*nc/4*(2*rand(ns,1)-1), pc*nr/4*(2*rand(ns,1)-1)]);
  S = S + r;
end

