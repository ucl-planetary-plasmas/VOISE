function a = homogeneousCriteria(W,varargin)
% function a = homogeneousCriteria(W,varargin)

%
% $Id: homogeneousCriteria.m,v 1.1 2009/02/08 21:07:18 patrick Exp $
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

if 0

  nW = length(W);

  if 0
    % mean and standard deviation
    m = mean(W);
    s = std(W);
  else
    % median and median absolute deviation
	  % are robust statistic estimators
	  % that are not unduly affected by small departures 
	  % from model assumptions
    m = median(W);
    s = mean(abs(W-m));
  end

  % find values outside range [m-s, m+s]
  ii = find(W < m-s | W > m+s);
  nOut = length(ii);

  % homogeneity criteria is proportion of pixels in W with value 
  % outside range [m-s, m+s];
  a = nOut/nW;

  if 0
  fprintf(1,'card(W in VR)=%d, card(W in VR & ~in [m-s,m+s])=%d, a=%2f\n', ...
	  nW, nOut, a);
  end

else

  % max-min criteria 
  D = (max(W)-min(W));
  % norm is max D(P), P in VD
  norm = varargin{1};

  a = D/norm;

end
