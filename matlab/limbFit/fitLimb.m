function fit = fitLimb(fit,Sx,Sy,Sls)
% function fit = fitLimb(fit,Sx,Sy,Sls)

%
% $Id: fitLimb.m,v 1.1 2009/10/16 14:02:32 patrick Exp $
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

% column vector of statistical weights 
if ~exist('Sls') | isempty(Sls),
  % constant
  W = ones(size(R));
else
  % proportional to 1/sqrt(var)
  W = 1./Sls(:);
end

% column vec of initial parameters
p0       = fit.p0;

% leasqr control parameters
stol     = fit.stol;
niter    = fit.niter;
dp       = fit.dp;
fracprec = fit.fracprec;
fracchg  = fit.fracchg;
options  = [fracprec, fracchg];

verbose  = fit.verbose;

[f,p,kvg,iter,corp,covp,covr,stdresid,Z,r2,ss] = ...
  leasqr(T, R, p0, 'limbModel', stol , niter, W, dp,'jacLimbModel',options);

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
if length(p) == 5,
  fprintf(1,'params est.: Xc=(%.1f,%.1f) a=%.1f b=%.1f tilt=%.0f\n', p);
  fprintf(1,'stdev  est.: Xc=(%.1g,%.1g) a=%.1g b=%.1g tilt=%.1g\n', psd)
elseif length(p) == 3,
  fprintf(1,'params est.: Xc=(%.1g,%.1g) R0=%.1f\n', p);
  fprintf(1,'stdev  est.: Xc=(%.1g,%.1g) R0=%.1f\n', psd)
end

