function checkVD(VD)
% function checkVD(VD)

%
% $Id: checkVD.m,v 1.1 2009/02/08 21:07:17 patrick Exp $
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

ok = true;
k = VD.k;

% check integrity of connectivity graph
for i = VD.Sk', % for all seeds at time k

  for j = VD.Nk{i}', % for all neighbours of seed with index i

    if isempty(find(VD.Nk{j} == i)), % 
		  s = sprintf('*** Connectivity error detected:\n');
		  s = [s sprintf('N%d(%d) = ', k, i)];
			s1 = printSet([], VD.Nk{i}); s = [s s1];
			s = [s sprintf('\n')];
		  s = [s sprintf('but %d ~in N%d(%d) = ', i, k, j)];
			s1 = printSet([],VD.Nk{j}); s = [s s1];
			s = [s sprintf('\n')];
			if 1
			  error(s);
			else
			  fprintf(1,'%s\n', s);
			end
			ok = false;
		end

	end

end

if ok,
	fprintf(1,'Voronoi Diagram connectivity ok\n');
end

% check integrity of referenced seeds
if ~isempty(setdiff(VD.Sk(:), unique(VD.Vk.lambda(:))))
  s = sprintf('*** Seed error detected:\n');
	s = [s sprintf('S(%d) = ', k)];
	s = [s sprintf('%d ', VD.Sk(:))];
	s = [s sprintf('\n')];
	s = [s sprintf('lambda = ')];
	s = [s sprintf('%d ', unique(VD.Vk.lambda(:)))];
	error(s);
end

