function testEllipseFit(ns,pc,p,p0)
% function testEllipseFit(ns,pc,p,p0)
%
% Try for example
% testEllipseFit
% testEllipseFit(500,0.2,[3.5,-8.5,340,250,30],[3,-8,340,250,60])
%

%
% $Id: testEllipseFit.m,v 1.1 2009/10/13 15:25:48 patrick Exp $
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

pstate = pause('query');
pause('off')

if nargin==0,
  testDriver(500,0.1,[3.5,-8.5,340,250,0],[-10,10,285,295,0]);
else
  testDriver(ns,pc,p,p0);
end


function testDriver(ns,pc,p,p0)

global verbose
verbose=[0 1 0];

Xc = p(1);
Yc = p(2);
a  = p(3);
b  = p(4);
t0 = p(5);

% init seed of Mersenne-Twister RNG
rand('twister',10);

t = 360*rand(1,ns);

Sx0 = Xc + a*cosd(t)*cosd(t0) - b*sind(t)*sind(t0);
Sy0 = Yc + a*cosd(t)*sind(t0) + b*sind(t)*cosd(t0);

Sx = Xc + a*(1+pc*(0.5-rand(1,ns))).*cosd(t)*cosd(t0)-...
          b*(1+pc*(0.5-rand(1,ns))).*sind(t)*sind(t0);
Sy = Yc + a*(1+pc*(0.5-rand(1,ns))).*cosd(t)*sind(t0)+...
          b*(1+pc*(0.5-rand(1,ns))).*sind(t)*cosd(t0);

%plot(Sx0,Sy0,'x',Sx,Sy,'o');

LSS = 2*sqrt((Sx0-Sx).^2+(Sy0-Sy).^2);
fprintf(1,'LSS min %f max %f\n', [min(LSS), max(LSS)]);

% fit all parameters
dp = ones(1,5);
fp = fitEllipse([],[],LSS,Sx,Sy,[1:length(LSS)],p0,dp);

fprintf(1,'exact  Xc(%.1f,%.1f) a=%.1f b=%.1f inclination=%.0f\n', p([1:5]));
fprintf(1,'guess  Xc(%.1f,%.1f) a=%.1f b=%.1f inclination=%.0f\n', p0([1:5]));
fprintf(1,'fitted Xc(%.1f,%.1f) a=%.1f b=%.1f inclination=%.0f\n', fp([1:5]));

hold on
plot(Sx0,Sy0,'ok');
hold off
h = get(gca,'children');
legend(h([1 2 4 3]),'data','data+noise','initial','fitted')


