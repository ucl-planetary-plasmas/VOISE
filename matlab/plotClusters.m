function params = plotClusters(CVD, params)
% function params = plotClusters(CVD, params)
%
% Plot clustering results from python code knn.py written
% by Jack Scantlebury for his MSc project at UCL in 2018, and adapted from
% the kNN-enhance algorithm (k Nearest Neighbour) and described in the paper 
% from Caiyan Jia et al.,
% Node Attribute-enhanced Community Detection in Complex Networks, 2017,
% Nature, Scientific Reports 7, doi: 10.1038/s41598-017-02751-8
% https://www.nature.com/articles/s41598-017-02751-8
%
% Example: 
%
%     % generate the VOISE segmentation with output in directory
%     % [voise.root '/share/output/north_proj]
%     webVOISEdemo1 
%
%     % run the python clustering algorithm knn.py 
%     global voise
%     command = [voise.root '/python/knn.py '];
%     datadir = [voise.root '/share/output/north_proj'];
%     [status,result] = system([command datadir]);
%     if status == 0,
%       disp(result);
%     end
% 
%     load([datadir '/voise.mat'])
%     plotClusters(CVD, params)

%
% $Id: plotClusters.m,v 1.2 2020/05/03 18:28:57 patrick Exp $
%
% Copyright (c) 2020 Patrick Guio <patrick.guio@gmail.com>
% All Rights Reserved.
%
% This program is free software; you can redistribute it and/or modify it
% under the terms of the GNU General Public License as published by the
% Free Software Foundation; either version 3 of the License, or (at your
% option) any later version.
%
% This program is distributed in the hope that it will be useful, but
% WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General
% Public License for more details.
%
% You should have received a copy of the GNU General Public License
% along with this program. If not, see <http://www.gnu.org/licenses/>.

clf
load(strcat(params.oDir, '/clusters.txt'))

axis equal;

pbaspect([1 1 1]);

activeX = CVD.Sx(CVD.Sk);
activeY = CVD.Sy(CVD.Sk);

xlim([0, CVD.nc])
ylim([0, CVD.nr])

% create image of clustered Voronoi regions (VR)
Wop = sopToWop(CVD, clusters);

% create image of length-scale within a VR defined as 
% the square root of the area of the polygon
[WLS, SLS] = getVDOp(CVD, params.W, 'sqrtLen');

cluster_count = max(clusters(:)) + 1;

colormap(spring(cluster_count));

x = params.x;
y = params.y;

% image of clustered VR
imagesc(x, y, Wop);
set(gca, 'YDir', 'normal')

hold on

W = CVD.W;
extentX = W.xM - W.xm;
extentY = W.yM - W.ym;

sx = (max(params.x) - min(params.x)) / extentX;
sy = (max(params.y) - min(params.y)) / extentY;

Sx = CVD.Sx(CVD.Sk);
Sy = CVD.Sy(CVD.Sk);

% calculate Voronoi diagram for seeds with coordinates (Sx(i),Sy(i))
[vx, vy] = voronoi(Sx, Sy);
Vx = (vx-W.xm)*sx+min(params.x);
Vy = (vy-W.ym)*sy+min(params.y);

plot(Vx, Vy, '-k', 'LineWidth', 0.5)
hold off

axis square;

c = colorbar;
c.Ticks = [];
c.TickLabels = c.Ticks;

%title("Clustering: " + "$\bar{s}($" + num2str(cluster_count) +...
%    "$) = 0.40$", 'Interpreter' ,'latex');
title("Clustering: " +"$k = $" +num2str(cluster_count), ...
      'Interpreter', 'latex');

xlabel(sprintf('x [%s]', params.pixelUnit{1}))
ylabel(sprintf('y [%s]', params.pixelUnit{2}))

print([params.oDir, 'clusters'], '-depsc');

%dpi = strcat('-r', num2str(params.dpi));
%print([params.oDir, 'clusters'], '-dpdf', dpi);
%print([params.oDir, 'clusters'], '-depsc', '-opengl', dpi);

