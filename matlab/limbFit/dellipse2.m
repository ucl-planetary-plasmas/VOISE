function dr=dellipse2(xy,f,p,dp,func)
% function dr=dellipse2(xy,f,p,dp,func)

%
% $Id: dellipse2.m,v 1.1 2010/07/14 15:42:36 patrick Exp $
%
% Copyright (c) 2010
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

global verbose

if ~isempty(verbose) & length(verbose)>2 & verbose(3),
  fprintf(1,'calling dellipse2 (xc,yc,a,b,t0)=(%.1f,%.1f,%.1f,%.1f,%.0f)\n',...
          p(1:5));
end


xc = p(1); % x-coordinate of ellipse center 
yc = p(2); % y-coordinate of ellipse center
a  = p(3); % semi-major axis
b  = p(4); % semi-minor axis
t0 = p(5); % tilt angle of semi-major axis to x-axis [rad]
ti = p(6:end); % angles

ni = length(xy);
m = fix(ni/2);

xi = xy(1:m);
yi = xy(m+1:ni);


ci = cos(ti);
si = sin(ti);

Q = rot(t0);
Qp = rotprime(t0);

dgdti = zeros(ni,1);
dgdt0 = zeros(ni,1);
dgda = zeros(ni,1);
dgdb = zeros(ni,1);
dgdxc = zeros(ni,1);
dgdyc = zeros(ni,1);
for i=1:m,
  % d/dti
  dgdti(i+[0,m]) = Q*[-a*si(i);b*ci(i)];
	% d/dt0
	dgdt0(i+[0,m]) = Qp*[a*ci(i);b*si(i)];
	% d/da
	dgda(i+[0,m]) = Q*[ci(i);0];
	% d/db
	dgdb(i+[0,m]) = Q*[0;si(i)];
	% d/dxc
	dgdxc(i+[0,m]) = [1;0];
	% d/dyc
	dgdyc(i+[0,m]) = [0;1];
end

% drs is [drs/dxc, drs/dyc, drs/da, drs/db, drs/dt0, drs/dti....]
dr = [dgdxc, dgdyc, dgda, dgdb, dgdt0,[diag(dgdti(1:m));diag(dgdti(m+1:ni))]];


function Q = rot(alpha)

Q = [cos(alpha), -sin(alpha); ...
     sin(alpha), cos(alpha)];

function Qp = rotprime(alpha)

Qp = [-sin(alpha), -cos(alpha); ...
      cos(alpha), -sin(alpha)];