function S = fig4Seeds(nr,nc,ns,varargin)
% function S = fig4Seeds(nr,nc,ns,['pc',value])
%
% string 'pc' followed by a value is an optional argument.
% pc is a percentage to indicate the relative fluctuation introduced
% in the randomisation of the regular tesselation (default pc = 0.25)

%
% $Id: fig4Seeds.m,v 1.2 2009/11/10 14:53:06 patrick Exp $
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



ns = 9;

x = round(nc*[1/4; 1/2; 3/4]);
y = round(nr*[1/4; 1/2; 3/4]);

% initialise array S(ns,2) 
% seed s has coordinates (x,y) = S(s, 1:2) 
S = [[x(1:3); x([1 3]); x(1:3); x(2)], ...
     [y(3)*ones(3,1); y(2)*ones(2,1); y(1)*ones(3,1); y(2)]];

pc = getfield(parseArgs(struct('pc',0.25), varargin{:}),'pc');


if pc, % random fluctuation of 100*pc % of distance between seeds
  r = round([pc*nc/4*(2*rand(ns,1)-1), pc*nr/4*(2*rand(ns,1)-1)]);
  S = S + r;
end

