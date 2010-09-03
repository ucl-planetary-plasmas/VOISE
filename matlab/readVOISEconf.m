function params = readVOISEconf(conffile)

%
% $Id: readVOISEconf.m,v 1.1 2010/09/03 17:12:33 patrick Exp $
%
% Copyright (c) 2010
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

% load default VOISE parameters
params = getDefaultVOISEParams();


% read VOISE configuration file
% all fields initialised by getDefaultVOISEParams 
% can be redefined in the configuration file
% The syntax is
% field = value
fid = fopen(conffile);
C = textscan(fid, '%s=%s');
fclose(fid);
key = C{1};
val = C{2};
v = str2double(val); idx = ~isnan(v);
val(idx) = num2cell(v(idx));


for i=1:length(key),
  if isfield(params, key{i}),
    params.(key{i}) = val{i};
  end
end


