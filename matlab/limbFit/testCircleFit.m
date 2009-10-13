function testCircleFit(ns,pc,p,p0)
% function testCircleFit(ns,pc,p,p0)
% 
% Try for example
% testCircleFit
% testCircleFit(100,0.2,[3.5,-8.5,250],[-20,10,285])
%


%
% $Id: testCircleFit.m,v 1.1 2009/10/13 16:04:59 patrick Exp $
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
  testDriver(500,0.05,[3.5,-8.5,250],[-10,10,285]);
elseif nargin==4,
  testDriver(ns,pc,p,p0);
end

pause(pstate)

function testDriver(ns,pc,p,p0)

global verbose
verbose=[0 1 0];

Xc = p(1);
Yc = p(2);
r  = p(3);

% init seed of Mersenne-Twister RNG
rand('twister',10);

t = 360*rand(1,ns);

Sx0 = Xc + r*cosd(t);
Sy0 = Yc + r*sind(t);

Sx = Xc + r*(1+pc*(0.5-rand(1,ns))).*cosd(t);
Sy = Yc + r*(1+pc*(0.5-rand(1,ns))).*sind(t);

%plot(Sx0,Sy0,'x',Sx,Sy,'o');

fprintf(1,'seeds # %d,  pc %f\n', ns, pc);

LSS = 2*sqrt((Sx0-Sx).^2+(Sy0-Sy).^2);
fprintf(1,'LSS min %f max %f\n', [min(LSS), max(LSS)]);

fp = fitCircle([],[],LSS,Sx,Sy,[1:length(LSS)],p0);

fprintf(1,'exact  Xc(%.1f,%.1f) R=%.1f\n', p([1:3]));
fprintf(1,'guess  Xc(%.1f,%.1f) R=%.1f\n', p0([1:3]));
fprintf(1,'fitted Xc(%.1f,%.1f) R=%.1f\n', fp([1:3]));

hold on
plot(Sx0,Sy0,'ok');
hold off
h = get(gca,'children');
legend(h([1 2 4 3]),'data','data+noise','initial','fitted')


