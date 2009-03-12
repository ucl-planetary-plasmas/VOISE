function [hs,hl]=plotDT(x,y,scatterSpecs,lineSpecs)
% function [hs,hl]=plotDT(x,y)

%
% $Id: plotDT.m,v 1.1 2009/03/12 15:41:24 patrick Exp $
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

tri = delaunay(x, y);

if ~exist('scatterSpecs','var') | isempty(scatterSpecs),
  hs = scatter(x, y);
else
  hs = scatter(x, y, scatterSpecs{:});
end

hl = zeros(size(tri,1));
if ~exist('lineSpecs','var') | isempty(lineSpecs),
  for i=1:size(tri,1)
    index = [tri(i,:) tri(i,1)];
	  hl(i) = line(x(index), y(index)); %,'color','k');
  end
else
  for i=1:size(tri,1)
    index = [tri(i,:) tri(i,1)];
	  hl(i) = line(x(index), y(index),lineSpecs{:}); %,'color','k');
  end
end

