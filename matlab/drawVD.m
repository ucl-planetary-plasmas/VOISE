function drawVD(VD)
% function drawVD(VD)

% Algorithm described in 
% Discrete Voronoi Diagrams and the SKIZ Operator: A Dynamic Algorithm
% R. E. Sequeira and F. J. Preteux
% IEEE Transactions on pattern analysis and machine intelligence
% Vol. 18, No 10, October 1997

%
% $Id: drawVD.m,v 1.1 2009/02/08 21:07:15 patrick Exp $
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

[vx,vy]=voronoi(VD.Sx(VD.Sk), VD.Sy(VD.Sk));

subplot(211)
K = VD.Vk.v; K(VD.Vk.v==1) = -VD.Vk.lambda(VD.Vk.v==1);
imagesc(VD.Vk.lambda + K)
axis xy
%colorbar
hold on
for i = VD.Sk',
  plot(VD.Sx(i), VD.Sy(i), 'xk', 'MarkerSize',8)
	text(VD.Sx(i), VD.Sy(i), num2str(i), 'verticalalignment', 'bottom');
end
plot(vx,vy,'-k','LineWidth',0.5)
hold off
title(sprintf('k = %d card(S) = %d', VD.k, length(VD.Sk)))

subplot(212), 
mu  = (VD.x-VD.Sx(VD.Vk.lambda)).^2+(VD.y-VD.Sy(VD.Vk.lambda)).^2;
imagesc(mu); 
axis xy
%colorbar
hold on
for i = VD.Sk',
  plot(VD.Sx(i), VD.Sy(i), 'xk', 'MarkerSize',8)
	text(VD.Sx(i), VD.Sy(i), num2str(i), 'verticalalignment', 'bottom');
end
plot(vx,vy,'-k','LineWidth',0.5)
hold off

drawnow
%pause
