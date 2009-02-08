function createTestImage
% function createTestImage

%
% $Id: createTestImage.m,v 1.1 2009/02/08 21:07:18 patrick Exp $
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

x = linspace(-1,1,200);
y = linspace(-1,1,200);
[X,Y]=meshgrid(x,y);

Z = 10*zeros(200,200);

a = 3/4; b = 3/4;
Z(X.^2/a^2+Y.^2/b^2 <= 1) = 20;

a = 1/2; b = 2/3;
Z(X.^2/a^2+Y.^2/b^2 < 1) = 30;

a = 1/3; b = 1/4;
Z(X.^2/a^2+Y.^2/b^2 < 1) = 20;

randn('state', 10);

Z = Z + randn(200,200);

Z = fix(255*(Z-min(Z(:)))/(max(Z(:))-min(Z(:))))+1;

imagesc(Z)

save ../share/testImage x y Z
