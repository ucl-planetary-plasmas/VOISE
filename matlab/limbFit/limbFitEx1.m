function [varargout] = limbFitEx1(action)
% function varargout = limbFitEx1(action)
%
% where action is any of 'clean', 'VOISE' and 'limbFit'
%
% limbFitEx1('clean');
% [params,IVD,DVD,MVD,CVD] = limbFitEx1('VOISE')
% [params,fit1,fit2] = limbFitEx1('limbFit')

%
% $Id: limbFitEx1.m,v 1.1 2009/11/06 17:15:23 patrick Exp $
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

switch lower(action),

   case 'clean',

	   unix('rm -rf ../../share/alljup_UT1013');

   case 'voise',

     % initial number of seeds
     ns = 30;

     % load default parameters
     params = getDefaultVOISEParams;

     % set image filename
     params.iFile = '../../share/alljup_UT1013.fits';

     % set output parameters
     params.oDir = '../../share/alljup_UT1013/';
     params.oMatFile = 'voise';

     % load image from fits file
     params.Wo = fitsread(params.iFile);

     % set up filter parameters
     params.winSize = [11,11];
     params.filter = @median;
     params.W  = filterImage(params.Wo,params.winSize,params.filter);
     %[params.W,xbins,cdf] = histEq(params.W,256);

     % set axes limits
     params.Wlim = [min(params.W(:)) max(params.W(:))];

     params.x0 = (1+size(params.W,2))/2;
     params.x = [1:size(params.W,2)]-params.x0;
     params.xlim = [min(params.x) max(params.x)];

     params.y0 = (1+size(params.W,1))/2;
     params.y = [1:size(params.W,1)]-params.y0;
     params.ylim = [min(params.y) max(params.y)];

     if ~exist(params.oDir,'dir')
       unix(['mkdir ' params.oDir]);
     end

     % VOISE parameters
     params.dividePctile    = 98;
     params.d2Seeds         = 10;
     params.mergePctile     = 50;
     params.dmu             = 0.2;
     params.thresHoldLength = 0.3;
     params.regMaxIter      = 2;

     % report parameters
     params.divideExport = false;
     params.mergeExport = false;
     params.movDiag = false;
     params.movPos=[600 50 1300 1050]; % movie window size

     % init seed of Mersenne-Twister RNG
     rand('twister',10);

     close all

     t = cputime;
     [params,IVD,DVD,MVD,CVD] = VOISE(params, ns, @randomSeeds);
     t = cputime-t;
     fprintf(1,'Elapsed time %2d:%2d:%2d\n', floor(t/3600), ...
             floor(mod(t,3600)/60), floor(mod(mod(t,3660),60)));

		 varargout{1} = {params};
		 varargout{2} = {IVD};
		 varargout{3} = {DVD};
		 varargout{4} = {MVD};
		 varargout{5} = {CVD};

   case 'limbfit',

	   path = '../../share/alljup_UT1013/';
		 data = 'voise';
		 fprintf(1,'Opening %s%s\n',path,data);
		 eval(['load ' path data]);


		 VD = DVD;
		 vdtype = 'DVD';

     % initial guess
		 if 1
		   p0 = [0 0 300 280 0];
     else
		   p0 = getEllipseParams();
     end
		 fit1 = getDefaultFitParams(p0);
		 fit1.LSmax = 16;
		 fit1.Rmin = 0.9;
		 fit1.Rmax = 1.1;
		 fit1.VD = vdtype;
		 fit1.dp = [1 1 1 1 0];
		 fit1 = getLimb(CVD, params, fit1);

		 % try again with better selected seeds
		 fit2 = getDefaultFitParams(fit1.p);
		 fit2.LSmax = 12;
		 fit2.Rmin = 0.95;
		 fit2.Rmax = 1.05;
		 fit2.VD = vdtype;
		 fit2.dp = [1 1 1 1 0];
		 fit2 = getLimb(CVD, params, fit2);

		 orient landscape, set(gcf,'PaperOrientation','portrait');
		 exportfig(gcf,[path '/select.eps'],...
		               'color','cmyk','boundscode','mcode','LockAxes',0);

     subplot(111)
		 plotLimbFit(params,fit1,fit2);

		 orient landscape, set(gcf,'PaperOrientation','portrait');
		 exportfig(gcf,[path '/fit.eps'],...
		           'color','cmyk','boundscode','mcode','LockAxes',0);

		 varargout{1} = {params};
		 varargout{2} = {fit1};
		 varargout{3} = {fit2};

end
