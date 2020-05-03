function [A,D,L,l,d,s,I] = VD2GraphMatrices(VD)
% function [A,D,L,l,d,s,I] = VD2GraphMatrices(VD)
%
% Create matrice representation of graph in VD
%
% A: Adjacency matrix
% D: Degree matrix
% L: Laplacian matrix
% l: Normalised Laplacian matrix
% d: Distance matrix
% s: Similarity matrix
% I: Incidence matrix
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
%     [A,D,L,l,d,s,I] = VD2GraphMatrices(CVD);
%
%     close all
%
%     figure
%     subplot(221) 
%     imagesc(A)
%     title('adjacency matrix')
%     colormap(gca,parula(2))
%     colorbar
%
%     subplot(222)
%     imagesc(D)
%     title('degree matrix')
%     colormap(gca,parula(max(D(:))-min(D(:))))
%     colorbar
%
%     subplot(223)
%     imagesc(L)
%     title('Laplacian matrix')
%     colormap(gca,parula(max(L(:))-min(L(:))))
%     colorbar
%
%     subplot(224)
%     imagesc(l)
%     title('Laplacian matrix')
%     colormap(gca,parula(256))
%     colorbar
%
%     figure
%     subplot(311)
%     imagesc(d)
%     title('distance matrix')
%     colormap(gca,parula(max(d(:))-min(d(:))))
%     colorbar
%
%     subplot(312)
%     imagesc(s)
%     title('similarity matrix')
%     colormap(gca,parula(256))
%     colorbar
%
%     subplot(313)
%     imagesc(I)
%     title('incidence matrix')
%     colormap(gca,parula(2))
%     colorbar


%
% $Id: VD2GraphMatrices.m,v 1.2 2020/05/03 21:10:23 patrick Exp $
%
% Copyright (c) 2015 Patrick Guio <patrick.guio@gmail.com>
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

% Adjacency matrix
A = zeros(length(VD.Nk));
% Degree matrix
D = zeros(length(VD.Nk));

% distance matrix
d = zeros(length(VD.Nk));
% Similarity matrix
s = zeros(length(VD.Nk));

for i = 1:length(VD.Nk), % for all seeds that ever been registered
  for n = 1:length(VD.Nk{i}),
	  j = VD.Nk{i}(n);
		% Adjacency
		A(i,j) = 1;
		% distance
		d(i,j) = sqrt((VD.Sx(i)-VD.Sx(j)).^2+(VD.Sy(i)-VD.Sy(j)).^2);
		% similarity
		s(i,j) = exp(-d(i,j).^2/2);
	end
	% degree
  D(i,i) = length(VD.Nk{i});
end

% Laplacian matrix
L = D-A;

% Normalised Laplacian matrix
% l = eye(size(L))-D^(-1/2)*A*D^(-1/2)
l = D^(-1/2)*L*D^(-1/2);

% Incidence matrix
nvertices = length(VD.Nk);
Aupper = triu(A);
nedges = sum(Aupper(:));
I = zeros(nvertices,nedges);
[ii,jj] = find(A==1);
for i = 1:nedges,
  I(ii(i),i) = 1;
  I(jj(i),i) = 1;
end
