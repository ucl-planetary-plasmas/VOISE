function seedDist(VD,params)
% function seedDist(VD,params)

%
% $Id: seedDist.m,v 1.5 2010/04/14 07:24:19 patrick Exp $
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

sx = (max(params.x)-min(params.x))/(VD.xM-VD.xm);
sy = (max(params.y)-min(params.y))/(VD.yM-VD.ym);

X = (VD.Sx(VD.Sk)-VD.xm)*sx+min(params.x);
Y = (VD.Sy(VD.Sk)-VD.ym)*sy+min(params.y);

subplot(111);
plot(X,Y,'o')
xlabel('x [R_J]')
ylabel('y [R_J]')
title('Seeds spatial distribution')

printFigure(gcf,[params.oDir 'seeddist1.eps']);

R = sqrt(X.^2+Y.^2);
T = atan2(Y,X)*180/pi;

clf 
subplot(111)
plot(R,T,'o')
h=ylabel('\theta [deg]','VerticalAlignment','top');
%set(h)
xlabel('\rho [R_J]')
title('Seeds spatial distribution')


ri = linspace(0,1.1*max(R),30);
ti = linspace(-90,90,80);

[Ri,Ti] = meshgrid(ri,ti);

Hi = zeros(size(Ri));

for i=1:length(ti)-1,
  for j=1:length(ri)-1,
	  Hi(i,j) = length(find(R>=Ri(i,j) & R<Ri(i+1,j+1) & T>=Ti(i,j) & T<Ti(i+1,j+1)));
	end
end

printFigure(gcf,[params.oDir 'seeddist2.eps']);

return

printFigure(gcf,[params.oDir 'seeddist.eps']);

return 

subplot(313)
%imagesc(ti,ri,Hi')
pcolor(Ti',Ri',Hi'), shading interp
axis xy
%colorbar
return

hist(T,ti),pause
hist(R,ri),pause
