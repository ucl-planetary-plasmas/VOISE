function printVOISEsetup(params)
% function printVOISEsetup(params)

%
% $Id: printVOISEsetup.m,v 1.2 2009/11/11 13:38:39 patrick Exp $
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

global voise

header = sprintf('VOISE Set Up -- version %s', voise.version);
line = sprintf('%s', char('*'*ones(length(header),1)));


fprintf(1,'\n%s\n%s\n%s\n\n', line, header, line);

fprintf(1,' * Dividing parameters\n');
fprintf(1,'   -------------------\n\n');
fprintf(1,'   p_D             = %.1f\n', params.dividePctile);
fprintf(1,'   d^2_m           = %.1f\n', params.d2Seeds);
fprintf(1,'   algo            = %d\n\n', params.divideAlgo);

fprintf(1,' * Merging parameters\n');
fprintf(1,'   -------------------\n\n');
fprintf(1,'   p_M             = %.1f\n', params.mergePctile);
fprintf(1,'   dmu             = %.1f\n', params.dmu);
fprintf(1,'   thresHoldLength = %.1f\n', params.thresHoldLength);
fprintf(1,'   algo            = %d\n\n', params.mergeAlgo);

fprintf(1,' * Regularising parameters\n');
fprintf(1,'   -----------------------\n\n');
fprintf(1,'   regMaxIter      = %d\n'  , params.regMaxIter);
fprintf(1,'   algo            = %d\n\n', params.regAlgo);

fprintf(1,'\n%s\n', line);


