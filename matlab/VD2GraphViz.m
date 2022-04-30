function varargout = VD2GraphViz(filename, VD, params)
% function [gvString]=VD2GraphViz(filename, VD, params)
%
% Create GraphViz graph file (.gv) from a Voronoi diagram and 
% returns on demand a string containing the graph file.
%
% Example:
%
%     % generate the VOISE segmentation with output in directory
%     % [voise.root '/share/output/north_proj]
%     webVOISEdemo1
%
%     global voise
%     datadir = [voise.root '/share/output/north_proj'];
%     load([datadir '/voise.mat'])
%     VD2GraphViz('CVD.gv',CVD, params);
%     

%
% $Id: VD2GraphViz.m,v 1.3 2022/04/30 15:45:25 patrick Exp $
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

% A list of the attributes available for GraphViz can be found at
% http://www.graphviz.org/content/attrs

% Put the GraphViz file into the output directory
filename = [params.oDir filename];

% Calculate equivalent scale length from VD as sqrt(S)
[imls, Sls] = getVDOp(VD, params.W, 'sqrtLen');

% Calculate median intensity from VD 
[immu, Smu] = getVDOp(VD, params.W, 'median');

% current time
k = VD.k;

s = sprintf('/*\n%s\n*/\n', ['sfdp ' ...
             '-Gsize="7,8" -Nshape=box -Efontsize=8 -Goverlap=prism ' ...
						 filename ' -Tpdf > graph.pdf']);

% If the graph is strict then multiple edges are not allowed between the
% same pairs of nodes
s = [s sprintf('\nstrict graph {\n')];

%s = [s sprintf('%s\n%s\n%s\n', s2, s1, s2)];

s1 = sprintf('\tgraph [splines=true overlap=false dim=3]');
s2 = sprintf('\trankdir=LR; // Left to Right, instead of Top to Bottom');

s = [s sprintf('%s\n%s\n', s1,s2)];

% Attributes Graphviz documentation 
% https://graphviz.org/doc/info/attrs.html

%for i = VD.Sk', % for all seeds at current time 
for i = 1:length(VD.Nk), % for all seeds that ever been registered
  s = [s sprintf('\t%d ', i)];
  s = [s sprintf('[pos="%.1f,%.1f,%.1f!" width=%.1f];\n',...
                 VD.Sx(i),VD.Sy(i),Smu(i),Sls(i))];
end
for i = 1:length(VD.Nk), % for all seeds that ever been registered
  s = [s sprintf('\t%d -- ', i)];
	sets = VD.Nk{i};
	s1 = printSet([], sets(sets>i)); s = [s s1];
  s = [s sprintf(';\n')];
end

s = [s sprintf('%s\n\n','}')];

fid = fopen(filename,'w');

if ~isempty(fid),
	fprintf(fid, '%s', s);
end

fclose(fid);

if nargout>0,
	varargout(1) = s;
end

