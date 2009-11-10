function s = parseArgs(s,varargin)
% function s = parseArgs(s,varargin)
%
% function to parse a list of arguments. If an element of the list is
% recognised as a field of the structure s, the next element of the list 
% is assigned to the field of the structure.

%
% $Id: parseArgs.m,v 1.1 2009/11/10 14:50:56 patrick Exp $
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
    s = setfield(s, varargin{i}, varargin{i+1});
	end
end