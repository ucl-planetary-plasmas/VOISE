function printFigure(hf,filename)
% function printFigure(hf,filename)

%
% $Id: printFigure.m,v 1.2 2009/11/12 17:33:04 patrick Exp $
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

if exist('exportfig','file') == 2,
  exportfig(hf, filename, 'color', 'cmyk', 'boundscode','mcode','LockAxes',0);
else
  print(hf, '-depsc', filename);
end

