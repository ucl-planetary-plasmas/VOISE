function s = parseArgs(s,varargin)
% function s = parseArgs(s,varargin)
%
% function to parse a list of arguments. If an element of the list is
% recognised as a field of the structure s, the next element of the list 
% is assigned to the field of the structure.

%
% $Id: parseArgs.m,v 1.3 2011/02/16 12:46:34 patrick Exp $
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

if isempty(varargin)
  return
end

for i=1:length(varargin)-1,
  if ischar(varargin{i}) & isfield(s,varargin{i}),
	  oldvalue = getfield(s, varargin{i});
		newvalue = varargin{i+1};
	  if isa(newvalue,'numeric') & isa(oldvalue,'numeric'),
		  if all(size(newvalue) == size(oldvalue)) | isempty(oldvalue), 
			% numeric arrays needs to have same size or variable originally empty
        s = setfield(s, varargin{i}, newvalue);
	    else
		    fprintf(1,'oldvalue: '); disp(oldvalue)
			  fprintf(1,'newvalue: '); disp(newvalue)
		    error(['New value for ''' varargin{i} ''' not adequate'])
		  end
		else
		  s = setfield(s, varargin{i}, newvalue);
		end
	end
end
