function varargout = getFitsPrimaryDataKeywordsValue(filename,keywords)
% function values = getFitsPrimaryDataKeywordsValue(filename,keywords)
%
% Return values of keywords contained in the Primary Data structure of a 
% file in FITS format.
% 
% keywords is a cell of strings representing the keywords such as
% {'PCX, 'PCY','CML','UDATE'}
% 
% values is a list of variable (one for each keyword). For the example
% given above values should look like
% [pcx,pcy,cml,udate]

%
% $Id: getFitsPrimaryDataKeywordsValue.m,v 1.2 2011/03/26 17:16:55 patrick Exp $
%
% Copyright (c) 2011 Patrick Guio <patrick.guio@gmail.com>
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

info=fitsinfo(filename);

values = cell(size(keywords));

for k = 1:length(keywords),

  keyword = keywords{k};
  for i = 1:size(info.PrimaryData.Keywords,1),

    if strcmp(info.PrimaryData.Keywords{i,1},keyword),

		  values{k} = info.PrimaryData.Keywords{i,2};

		end

	end

end

varargout = {values{:}};
