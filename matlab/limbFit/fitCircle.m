function p = fitCircle(VD,params,LSS,Sx,Sy,ii,p0)
% function p = fitCircle(VD,params,LSS,Sx,Sy,ii,p0)

%
% $Id: fitCircle.m,v 1.2 2010/01/12 12:40:55 patrick Exp $
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

global datapath ploc

% convert cartesian into polar coordinates
R = sqrt(Sx(ii).^2+Sy(ii).^2)';
T = 180/pi*atan2(Sy(ii),Sx(ii))';

fracprec=[0; 0; 0;];
fracchg=[Inf; Inf; Inf];
options=[fracprec fracchg];

stol = 1e-8;
niter=50;
dp=[1 1 1];

W = ones(size(R));
W = (1./LSS(ii)');

[f,p,kvg,iter,corp,covp,covr,stdresid,Z,r2,ss]=...
  leasqr(T, R, p0, 'circle', stol , niter, W, dp,...
	  'dcircle',options);

% degrees of freedom
nu = length(f) - length(p(dp==1));
% reduced chi2 statistic
chi2 = ss/(nu-1);
% standard deviation for estimated parameters 
psd = zeros(size(p));
sd = sqrt(diag(covp));
psd(dp==1) = sd;

fprintf(1,'conv %d iter %d r2 %.2f chi2 %f\n', kvg, iter, r2, chi2);
fprintf(1,'params est.: Xc(%.1f,%.1f) R=%.1f\n', p);
fprintf(1,'stdev  est.: Xc(%.1f,%.1f) R=%.1f\n', psd);

td = linspace(-180,180,50);

r0 = circle(td,p0);
r = circle(td,p);

%subplot(212)
clf
if ~isempty(params),
  figure(1)
	clf
else
  subplot(211)
end

plot(T,R,'o',td,r0,'-',td,r,'-');
set(gca,'ylim',[0 350],'xlim',[-180 180]);
xlabel('\theta [deg]');
ylabel('\rho [pixels]')

[legh,objh,oh,om] = legend('data','initial','fitted','location','SouthEast');
set(objh(1),'fontsize',9);

if ~isempty(datapath) & ~isempty(ploc),

  orient landscape, set(gcf,'PaperOrientation','portrait');
  exportfig(gcf,[datapath '/phasloc' num2str(ploc) '.eps'],...
            'color','cmyk','boundscode','mcode','LockAxes',0);

  ploc = ploc+1;

end


%subplot(222)

if ~isempty(params),
  imagesc(params.x,params.y,params.W);
  axis xy
  axis equal
  axis tight
  xlabel('x [pixels]')
  ylabel('y [pixels]')
	hold on
else
  figure(2)
	clf
	hold on
end

Ts = sort(T);
r0 = circle(Ts,p0);
r = circle(Ts,p);

Xe = r.*cosd(Ts);
Ye = r.*sind(Ts);

Xe0 = r0.*cosd(Ts);
Ye0 = r0.*sind(Ts);


% initial guess and fitted 
plot(Xe0,Ye0,'r--o', Xe,Ye, 'r-o','linewidth',1)
% data points
plot(Sx(ii), Sy(ii), 'go');

if ~isempty(VD) & ~isempty(params),
  x = params.x;
	y = params.y;
  line([p(1), p(1)],[min(params.y),max(params.y)],'color','w');
  line([min(params.x),max(params.x)],[p(2), p(2)],'color','w');
  [vx,vy]=voronoi(VD.Sx(VD.Sk), VD.Sy(VD.Sk));
  vx = (vx-min(VD.Sx))/(max(VD.Sx)-min(VD.Sx))*(max(x)-min(x))+min(x);
  vy = (vy-min(VD.Sy))/(max(VD.Sy)-min(VD.Sy))*(max(y)-min(y))+min(y);
  plot(vx,vy,'-k','LineWidth',0.5)
end
title(sprintf('X_c(%.1f,%.1f) R=%.1f', p([1:3])))

hold off

if ~isempty(datapath) & ~isempty(ploc),

  orient landscape, set(gcf,'PaperOrientation','portrait');
  exportfig(gcf,[datapath '/phasloc' num2str(ploc) '.eps'],...
            'color','cmyk','boundscode','mcode','LockAxes',0);

  ploc = ploc+1;

end
