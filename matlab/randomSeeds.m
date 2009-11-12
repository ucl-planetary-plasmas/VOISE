function [S,varargout] = randomSeeds(nr,nc,ns,varargin)
% function [S[,pc]] = randomSeeds(nr,nc,ns[,'pc',pc])
% 
% string 'pc' followed by a value for pc is optional.
% pc is a percentage to indicate the relative fluctuation introduced
% in the randomisation of the regular tesselation (default pc = 0.02)

%
% $Id: randomSeeds.m,v 1.5 2009/11/12 15:10:29 patrick Exp $
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

pc = getfield(parseArgs(struct('pc',0.02), varargin{:}),'pc');

% initialise array S(ns,2) 
% seed s has coordinates (x,y) = S(s, 1:2) 

% no seeds in the 100*pc % from the boundary of the image
S = round([(nc-2*pc*nc)*(rand(ns,1))+pc*nc, ...
           (nr-2*pc*nr)*(rand(ns,1))+pc*nr]);

% iterate until all seeds are different
uniqueSeeds=false;
while ~uniqueSeeds,

  nIdentical = 0;
  for k=1:size(S,1),
    ii = find(S(k,1)==S([k+1:end-1],1) & S(k,2)==S([k+1:end-1],2));
    if ~isempty(ii),
      ns = length(ii);
      S(k+ii,:) = round([(nc-2*pc*nc)*(rand(ns,1))+pc*nc, ...
                         (nr-2*pc*nr)*(rand(ns,1))+pc*nr]);
      nIdentical = nIdentical+ns;
    end
  end
  if ~nIdentical,
    uniqueSeeds = true;
  else
    fprintf(1,'nIdentical %d\n', nIdentical);
  end

end

if nargout > 1,
  varargout{1} = pc;
end
