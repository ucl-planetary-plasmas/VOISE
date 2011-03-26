function movieHandler(params,action,varargin)
% function movieHandler(params,action[,options])

%
% $Id: movieHandler.m,v 1.2 2011/03/26 17:16:55 patrick Exp $
%
% Copyright (c) 2009-2011 Patrick Guio <patrick.guio@gmail.com>
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

persistent aviobj

switch lower(action)

  case 'init',
	  set(gcf,'position',params.movPos);
		set(gcf,'DoubleBuffer','on');
		aviobj = avifile([params.oDir params.oMovFile],'fps',2);
	
	case 'addframe',
	  aviobj = addframe(aviobj, getframe(gcf,[0 0 params.movPos(3:4)]));
	
	case 'close',
	  aviobj = close(aviobj);

end
