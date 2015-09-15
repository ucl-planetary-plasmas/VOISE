function [xll,yll,xld,yld,xtl,ytl,xtd,ytd,xc,yc]=ltc(r,e,obs,sun)
% function [xll,yll,xld,yld,xtl,ytl,xtd,ytd,xc,yc]=ltc(r,e,obs,sun)

%
% $Id: ltc.m,v 1.4 2015/09/15 18:01:11 patrick Exp $
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

theta = linspace(0,2*pi,5000);


% deg into rad
olat = obs.lat*pi/180;
olon = obs.lon*pi/180;
slat = sun.lat*pi/180;
slon = sun.lon*pi/180;

dlon = slon-olon;

% planet 
xp = r*cos(theta);
yp = r*sqrt(1-e^2)*sin(theta);

% limb
xl = r*cos(theta);
yl = r*sqrt(1-e^2*cos(olat)^2)*sin(theta);

% cos and sin of latitude of normal to limb
csnllat = (1-e^2)*cos(olat)/sqrt((1-e^2)^2*cos(olat)^2+sin(olat)^2);
snnllat = sin(olat)/sqrt((1-e^2)^2*cos(olat)^2+sin(olat)^2);

%180/pi*acos(csnllat)
%180/pi*asin(snnllat)

DL = cos(slat)*sin(dlon)*xl+...
    (-cos(slat)*cos(dlon)*sin(olat)+sin(slat)/(1-e^2)*cos(olat))*yl;


xll = xl(DL>0);
yll = yl(DL>0);
xld = xl(DL<0);
yld = yl(DL<0);

% terminator seen from sun
xt0 = r*cos(theta);
yt0 = r*sqrt(1-e^2*cos(slat)^2)*sin(theta);

csntlat = (1-e^2)*cos(slat)/sqrt((1-e^2)^2*cos(slat)^2+sin(slat)^2);
snntlat = sin(slat)/sqrt((1-e^2)^2*cos(slat)^2+sin(slat)^2);

a = cos(dlon);
b = -snntlat*sin(dlon);
c = sin(dlon)*sin(olat);
d = snntlat*cos(dlon)*sin(olat)+csntlat*cos(olat);

A = [a, b; c, d];

xis = [0 A(1,1)];
yis = [0 A(2,1)];

xjs = [0 A(1,2)];
yjs = [0 A(2,2)];

xt = A(1,1)*xt0+A(1,2)*yt0;
yt = A(2,1)*xt0+A(2,2)*yt0;

% cusp points 
u = (1-e^2)*cos(slat)*cos(dlon)*sin(olat)-sin(slat)*cos(olat);
v = (1-e^2)*cos(slat)*sin(dlon);

t1 = sqrt(r^2*(1-e^2*cos(olat)^2)/(u^2*(1-e^2*cos(olat)^2)+v^2));
xc = u*[t1;-t1];
yc = v*[t1;-t1];

eL = eL(e,slat);
t2 = sqrt(r^2*(1-eL^2)*(a*d-b*c)^2/((d*v+b*u)^2*(1-eL^2)+(c*v+a*u)^2));

x2 = -v*[t2;-t2];
y2 = u*[t2;-t2];

thetaT = -atan(v/u);
aT = sqrt(u^2+v^2)*t1;
bT = sqrt(u^2+v^2)*t2;
fprintf(1,'aT=%.2f bT=%.2f thetaT=%.0f\n', aT, bT, 180/pi*thetaT);

% rotation \theta_T
rot = [cos(thetaT), sin(thetaT); -sin(thetaT), cos(thetaT)];

t = theta;
xt = aT*cos(t)*rot(1,1)+bT*sin(t)*rot(1,2);
yt = aT*cos(t)*rot(2,1)+bT*sin(t)*rot(2,2);
DT = csntlat*sin(dlon)*xt+(-csntlat*cos(dlon)*sin(olat)+snntlat*cos(olat))*yt;

plot(theta,DL,theta,DT);
pause

t=pi/2;
xt = aT*cos(t)*rot(1,1)+bT*sin(t)*rot(1,2);
yt = aT*cos(t)*rot(2,1)+bT*sin(t)*rot(2,2);
DT = csntlat*sin(dlon)*xt+(-csntlat*cos(dlon)*sin(olat)+snntlat*cos(olat))*yt;

DLC1 = cos(slat)*sin(dlon)*xc(1)+...
    (-cos(slat)*cos(dlon)*sin(olat)+sin(slat)*cos(olat)/(1-e^2))*yc(1);
DLC2 = cos(slat)*sin(dlon)*xc(2)+...
    (-cos(slat)*cos(dlon)*sin(olat)+sin(slat)*cos(olat)/(1-e^2))*yc(2);

DTC1 = csntlat*sin(dlon)*xc(1)+...
    (-csntlat*cos(dlon)*sin(olat)+snntlat*cos(olat))*yc(1);
DTC2 = csntlat*sin(dlon)*xc(2)+...
    (-csntlat*cos(dlon)*sin(olat)+snntlat*cos(olat))*yc(2);

fprintf(1,'DLC1=%.2g DLC2=%.2g DTC1=%.2g DTC2=%.2g\n', DLC1, DLC2, DTC1, DTC2);

xtl = xt(DT>0);
ytl = yt(DT>0);
xtd = xt(DT<0);
ytd = yt(DT<0);

if 0
if DT>0,
t = linspace(0,pi,100);
xtl = aT*cos(t)*rot(1,1)+bT*sin(t)*rot(1,2);
ytl = aT*cos(t)*rot(2,1)+bT*sin(t)*rot(2,2);
t = linspace(-pi,0,100);
xtd = aT*cos(t)*rot(1,1)+bT*sin(t)*rot(1,2);
ytd = aT*cos(t)*rot(2,1)+bT*sin(t)*rot(2,2);
else
t = linspace(-pi,0,100);
xtl = aT*cos(t)*rot(1,1)+bT*sin(t)*rot(1,2);
ytl = aT*cos(t)*rot(2,1)+bT*sin(t)*rot(2,2);
t = linspace(0,pi,100);
xtd = aT*cos(t)*rot(1,1)+bT*sin(t)*rot(1,2);
ytd = aT*cos(t)*rot(2,1)+bT*sin(t)*rot(2,2);
end
end

plot(xc,yc,'o',xll,yll,'-',xld,yld,'--');%,xtl,ytl,'-',xtd,ytd,'--');
     %xtl,ytl,'-',xtd,ytd,'--',xc,yc,'o'); %,x2,y2,'-o',xis,yis,xjs,yjs);

if 0
plot(xl,yl,xtl,ytl,'-',...
     xtd,ytd,'--',xc,yc,'-o',x2,y2,'-o',xis,yis,xjs,yjs);

[legh,objh,outh,outm] = legend('limb','tl','td','c','c2','is','js');
set(objh(1),'fontsize',9);

axis equal

title(sprintf('\\beta_{obs}=%.1f \\beta_{sun}=%.1f \\Delta\\lambda=%.1f',...
      olat*180/pi, slat*180/pi, dlon*180/pi));
end

function e = eL(e,beta)

e = e*sqrt(1-sin(beta)^2/(1-e^2*cos(beta)^2));
