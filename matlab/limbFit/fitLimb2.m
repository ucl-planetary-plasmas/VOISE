function fit = fitLimb2(fit,Sx,Sy,Sls)
% function fit = fitLimb2(fit,Sx,Sy,Sls)

%
% $Id: fitLimb2.m,v 1.6 2010/10/04 17:12:26 patrick Exp $
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

global verbose

% seeds position from Cartesian coordinates into polar coordinates
% column vector of observed values
R = sqrt(Sx(:).^2+Sy(:).^2);
% column vector or matrix of independent variables
T = 180/pi*atan2(Sy(:),Sx(:));

% arrange as [X;Y]
XY = [Sx(:); Sy(:)];

m = length(Sx(:));

% column vector of statistical weights 
if ~exist('Sls') | isempty(Sls),
  % constant
  W = [ones(size(R));ones(size(R))];
else
  % proportional to 1/sqrt(var)
  W = sqrt(2)*[1./Sls(:);1./Sls(:)];
end

% column vec of initial parameters
p0       = [fit.p0(:); T(:)*pi/180];

if length(fit.p0)==3,
  [xc,yc,R,a] = circfit(Sx,Sy);
	fprintf(1,'* circfit            Xc(%8.1f,%8.1f) R=%8.1f\n', xc,yc,R);
	Par = CircleFitByTaubin([Sx(:),Sy(:)]);
	fprintf(1,'* CircleFitByTaubin  Xc(%8.1f,%8.1f) R=%8.1f\n', Par);
elseif length(fit.p0)==5,
  p = ellipse_fit(Sx, Sy);
	fprintf(1,'* ellipse_fit        Xc(%8.1f,%8.1f) a=%8.1f b=%8.1f tilt=%8.2f\n', p);
	[A,p] = EllipseFitByTaubin([Sx(:),Sy(:)]);
	fprintf(1,'* EllipseFitByTaubin Xc(%8.1f,%8.1f) a=%8.1f b=%8.1f tilt=%8.2f\n', p);
	[A,p] = EllipseDirectFit([Sx(:),Sy(:)]);
	fprintf(1,'* EllipseDirectFit   Xc(%8.1f,%8.1f) a=%8.1f b=%8.1f tilt=%8.2f\n', p);
end

% leasqr control parameters
stol     = fit.stol;
niter    = fit.niter;
dp       = [fit.dp(:); ones(m,1)] ;
fracprec = [fit.fracprec; zeros(m,1)];
fracchg  = [fit.fracchg; Inf*ones(m,1)];
options  = [fracprec, fracchg];

verbose  = fit.verbose;

if length(fit.p0)==3,
[f,p,kvg,iter,corp,covp,covr,stdresid,Z,r2,ss] = ...
  leasqr(XY, XY, p0, 'circle2', stol , niter, W, dp,'dcircle2',options);
elseif length(fit.p0)==5,
[f,p,kvg,iter,corp,covp,covr,stdresid,Z,r2,ss] = ...
  leasqr(XY, XY, p0, 'ellipse2', stol , niter, W, dp,'dellipse2',options);
end

% degrees of freedom
nu = length(f) - length(p(dp==1));
% reduced chi2 statistic
chi2 = ss/(nu-1);
% standard deviation for estimated parameters 
psd = zeros(size(p));
sd = sqrt(diag(covp));
psd(dp==1) = sd;

% embed data in fit structure
fit.t = T;
fit.r = R;
fit.w = W;

% embed fit results in fit structure
fit.status = kvg;      % = 1 if convergence, = 0 otherwise
fit.iter   = iter;     % number of iterations used
fit.p      = p;        % trial or final parameters. i.e, the solution
fit.psd    = psd;      % standard deviation for estimated parameters
fit.corp   = corp;     % correlation matrix for parameters
fit.covp   = covp;     % covariance matrix of the parameters
fit.covr   = covr;     % diag(covariance matrix of the residuals)
fit.stdres = stdresid; % standardized residuals (res-m)/sd
fit.Z      = Z;        % matrix that defines confidence region 
fit.r2     = r2;       % coefficient of multiple determination
fit.ss     = ss;       % scalar sum of squares=sum-over-i(wt(i)*(y(i)-f(i)))^2
fit.nu     = nu;       % degrees of freedom
fit.chi2   = chi2;     % reduced chi2 statistic

fprintf(1,'status %d iter %d r2 %.2f chi2 %f\n', kvg, iter, r2, chi2);
if length(fit.p0) == 5,
  fprintf(1,'params estimated     Xc(%8.1f,%8.1f) a=%8.1f b=%8.1f tilt=%8.4f\n', p(1:5));
  fprintf(1,'stdev  estimated     Xc(%8.1f,%8.1f) a=%8.1f b=%8.1f tilt=%8.4f\n', psd(1:5))
elseif length(fit.p0) == 3,
  fprintf(1,'params estimated     Xc(%8.1f,%8.1f) R=%8.1f\n', p(1:3));
  fprintf(1,'stdev  estimated     Xc(%8.1f,%8.1f) R=%8.1f\n', psd(1:3))
end

phis = p(length(fit.p0)+1:end);
fprintf(1,'                     phi mean, min, max     %8.4f %8.4f %8.4f\n', ...
        mean(phis), min(phis), max(phis) );
phis = psd(length(fit.p0)+1:end);
fprintf(1,'                    dphi mean, min, max     %8.4f %8.4f %8.4f\n', ...
        mean(phis), min(phis), max(phis) );


if length(fit.p0)==3,
  % compute r
  xy = circle2(XY,fit.p);
	x = xy(1:m)-fit.p(1); y = xy(m+1:m+m)-fit.p(2);
	r = sqrt(x.^2+y.^2);
	% compute distance to seed
	d2seed = sqrt((x-Sx).^2+(y-Sy).^2);
	% compute dr
	xy = circle2(XY,fit.p-[0;0;fit.psd(3);zeros(m,1)]);
	x = xy(1:m)-fit.p(1); y = xy(m+1:m+m)-fit.p(2);
	rm = sqrt(x.^2+y.^2);
	xy = circle2(XY,fit.p+[0;0;fit.psd(3);zeros(m,1)]);
	x = xy(1:m)-fit.p(1); y = xy(m+1:m+m)-fit.p(2);
	rM = sqrt(x.^2+y.^2);
	dr = rM-rm;
	% compute dtheta
	xy = circle2(XY,fit.p-[zeros(3,1);fit.psd(4:end)]);
	x = xy(1:m)-fit.p(1); y = xy(m+1:m+m)-fit.p(2);
	tm = atan2(y,x);
	xy = circle2(XY,fit.p+[zeros(3,1);fit.psd(4:end)]);
	x = xy(1:m)-fit.p(1); y = xy(m+1:m+m)-fit.p(2);
	tM = atan2(y,x);
	dt = tM-tm;
	% fix the angles that might be wrongly interpreted
	dt(dt<0) = dt(dt<0)+2*pi;
	% compute infinitesimal volume rdrdtheta
  dS = r.*dr.*dt;
elseif length(fit.p0)==5,
  % compute r
  xy = ellipse2(XY,fit.p);
	x = xy(1:m)-fit.p(1); y = xy(m+1:m+m)-fit.p(2);
	r = sqrt(x.^2+y.^2);
	% compute distance to seed
	d2seed = sqrt((x-Sx).^2+(y-Sy).^2);
	% compute dr
	xy = ellipse2(XY,fit.p-[0;0;fit.psd(3:4);0;zeros(m,1)]);
	x = xy(1:m)-fit.p(1); y = xy(m+1:m+m)-fit.p(2);
	rm = sqrt(x.^2+y.^2);
	xy = ellipse2(XY,fit.p+[0;0;fit.psd(3:4);0;zeros(m,1)]);
	x = xy(1:m)-fit.p(1); y = xy(m+1:m+m)-fit.p(2);
	rM = sqrt(x.^2+y.^2);
	dr = rM-rm;
	% compute dtheta
	xy = ellipse2(XY,fit.p-[zeros(5,1);fit.psd(6:end)]);
	x = xy(1:m)-fit.p(1); y = xy(m+1:m+m)-fit.p(2);
	tm = atan2(y,x);
	xy = ellipse2(XY,fit.p+[zeros(5,1);fit.psd(6:end)]);
	x = xy(1:m)-fit.p(1); y = xy(m+1:m+m)-fit.p(2);
	tM = atan2(y,x);
	dt = tM-tm;
	% fix the angles that might be wrongly interpreted
	dt(dt<0) = dt(dt<0)+2*pi;
	% compute infinitesimal volume rdrdtheta
  dS = r.*dr.*dt;
end

figure
plot(Sls, sqrt(dS),'o',Sls,d2seed,'s')
axis equal
xlabel('Sls')
legend('sqrt(dS)','d2seed')
pause
close
