function [S,VDlim] = randomSeeds(nr,nc,ns,clipping)
% function [S,VDlim] = randomSeeds(nr,nc,ns,clipping)
%
% clipping is defined as percentage of the image size from each edge, i.e.
% the vector of length four with [left,right,bottom,top]
% default is a 2% default clipping from all edge [left,right,bottom,top]

%
% $Id: randomSeeds.m,v 1.8 2015/02/11 15:43:26 patrick Exp $
%
% Copyright (c) 2008-2012 Patrick Guio <patrick.guio@gmail.com>
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
% You should have received a copy of the GNU General Public License
% along with this program. If not, see <http://www.gnu.org/licenses/>.

if ~exist('clipping','var') || isempty(clipping),
  % 2% default clipping from all edge [left,right,bottom,top]
  clipping = [2, 2, 2, 2];
end
pc = clipping/100;

% initialise array S(ns,2) 
% seed s has coordinates (x,y) = S(s, 1:2) 
% where 1 < x < nc and 1 < y < nr
% i.e. no seeds on the border of the image

xm = floor(2 + (nc-3) * pc(1));
xM = ceil(2 + (nc-3) * (1-pc(2)));
ym = floor(2 + (nr-3) * pc(3));
yM = ceil(2 + (nr-3) * (1-pc(4)));

%fprintf(1,'(xm, xM, ym, yM) = (%d, %d, %d, %d)\n',xm,xM,ym,yM);

% uniform distribution over open range (boundary values not included)
%S = round([randraw('uniform', [xm, xM], ns, 1), ...
%           randraw('uniform', [ym, yM], ns, 1)]);
S = round([xm+(xM-xm)*rand(ns, 1), ...
           ym+(yM-ym)*rand(ns, 1)]);

% iterate until all seeds are different
uniqueSeeds=false;
while ~uniqueSeeds,

  nIdentical = 0;
  for k=1:size(S,1),
    ii = find(S(k,1)==S([k+1:end-1],1) & S(k,2)==S([k+1:end-1],2));
    if ~isempty(ii),
      ns = length(ii);
%      S(k+ii,:) = round([randraw('uniform', [xm, xM], ns, 1), ...
%                         randraw('uniform', [ym, yM], ns, 1)]);
      S(k+ii,:) = round([xm+(xM-xm)*rand(ns, 1), ...
                         ym+(yM-ym)*rand(ns, 1)]);
      nIdentical = nIdentical+ns;
    end
  end
  if ~nIdentical,
    uniqueSeeds = true;
  else
    fprintf(1,'nIdentical %d\n', nIdentical);
  end

end

% initialise VD seed limit structure
VDlim.xm = xm;
VDlim.xM = xM;
VDlim.ym = ym;
VDlim.yM = yM;
