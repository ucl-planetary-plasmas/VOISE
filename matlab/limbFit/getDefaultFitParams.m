function fit = getDefaultFitParams(p0,polePos)
% function fit = getDefaultFitParams(p0,polePos)

%
% $Id: getDefaultFitParams.m,v 1.3 2009/12/02 21:37:19 patrick Exp $
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

% initial parameters
fit.p0             = p0(:);
% default selection function for angles
fit.selectAngles   = @allAngles;
% allow to provide position of one of the poles
if nargin>1,
  fit.polePos      = polePos;
else
  fit.polePos      = [];
end

% leasqr control parameters

% scalar tolerance on fractional improvement in
% scalar sum of squares = sum((wt.*(y-f))^2); 
fit.stol           = 1e-6;

% scalar maximum number of iterations
fit.niter          = 50;

% desired fractional precision in parameter estimates
fit.fracprec       = zeros(length(p0), 1);

% maximum fractional step change in parameter vector
fit.fracchg        = Inf*ones(length(p0), 1);

% fractional increment of p for numerical partial derivatives
% i.e. when Jacobian not explicitly provided 
% dp(i) > 0 means central differences
% dp(i) < 0 means one-sided differences
% dp(i) 0 holds p(j) fixed
% Here jacobian of function 'limbModel' is provided by 'jacLimbModel' so that
% dp(i) != 0 means parameter fitted
fit.dp             =  ones(length(p0), 1);

% fit verbosity
% verbose(1): output progress report from leasqr (0: no, 1: yes)
% verbose(2): plot   progress report from leasqr (0: no, 1: yes)
% verbose(3): output from limb and derivative functions  (0: no, 1: yes)
fit.verbose        = [0 0 0];

% largest seed scale length to be considered for the fit of the limb
% *must* be initialised but no default value
fit.LSmax          = [];

% These values are used to select seeds for the fit of the limb
% the limb (as given by an ellipse/circle) is written in
% canonical form  
% x^2/a^2 + y^2/b^2 = 1
% The selected seeds must lie within the corona defined by 
% x^2/a^2 + y^2/b^2 >= Rmin^2 and x^2/a^2 + y^2/b^2 <= Rmax^2
fit.Rmin           = 0.9;
fit.Rmax           = 1.1;

