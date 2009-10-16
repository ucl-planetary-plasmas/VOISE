function fit = testLimbFit(ns,pc,p,p0,dp)
% function fit = testLimbFit(ns,pc,p,[p0],[dp])
% 
%    ns : number of seeds
%    pc : noise parameter
%    p  : circle/ellipse parameters
%    p0 : initial guess (optional)
%    dp : parameter fit flag  (0 means fixed, 1 means fitted)
%    
% Examples:
%
%  circle model
%  ------------
%
%  fit = testLimbFit(500,0.05,[3.5,-8.5,250],[-10,10,285]);
%  fit = testLimbFit(100,0.2,[3.5,-8.5,250],[-20,10,285]);
%  fit = testLimbFit(250,0.25,[-10.5,8.5,250]);
% 
%  ellipse model
%  -------------
%
%  fit = testLimbFit(500,0.1,[3.5,-8.5,340,250, 0],[-10,10,285,295,0 ]);
%  fit = testLimbFit(250,0.2,[3.5,-8.5,340,250,30],[  3,-8,340,250,60]);
%  fit = testLimbFit(250,0.2,[3.5,-8.5,340,250,30],...
%                    [  3,-8,340,250,30],[1 1 1 1 0]);
%  fit = testLimbFit(150,0.2,[3.5,-8.5,340,250,30]);

%
% $Id: testLimbFit.m,v 1.1 2009/10/16 14:13:44 patrick Exp $
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

pstate = pause('query');
pause('off')

if ~exist('p0'), 
  p0 = []; 
end

if ~exist('dp'), 
  dp = []; 
end

fit = testDriver(ns,pc,p,p0,dp);

pause(pstate)

function fit = testDriver(ns,pc,p,p0,dp)

% init seed of Mersenne-Twister RNG
rand('twister',10);

% generate random angles
t = 360*rand(1,ns);

% generate synthetic data
[Sx0,Sy0,Sx,Sy] = createLimbData(t,ns,p,pc);

if isempty(p0),
  clf
  plot(Sx,Sy,'o')
	if length(p) == 3,
    p0 = getCircleParams();
	elseif length(p) == 5,
    p0 = getEllipseParams();
	end
end

% mean average distance between exact and noisy data
LS = sqrt(((Sx0-Sx).^2+(Sy0-Sy).^2)/2);
fprintf(1,'LSS min %f max %f\n', [min(LS), max(LS)]);

fit = getDefaultFitParams(p0);
fit.verbose=[0 1 0];

if ~isempty(dp),
  if length(dp) == length(p0),
    fit.dp = dp;
	else
	  error('dp should be the same size as p0');
	end
end 


fit = fitLimb(fit,Sx,Sy,LS);

% print results
if length(p0) == 3,
  fprintf(1,'exact  Xc=(%.1f,%.1f) R=%.1f\n', p);
  fprintf(1,'guess  Xc=(%.1f,%.1f) R=%.1f\n', fit.p0);
  fprintf(1,'fitted Xc=(%.1f,%.1f) R=%.1f\n', fit.p);
else
fprintf(1,'exact  Xc=(%.1f,%.1f) a=%.1f b=%.1f tilt=%.0f\n', p);
fprintf(1,'guess  Xc=(%.1f,%.1f) a=%.1f b=%.1f tilt=%.0f\n', fit.p0);
fprintf(1,'fitted Xc=(%.1f,%.1f) a=%.1f b=%.1f tilt=%.0f\n', fit.p);
end

% sorted angles theta for data
ts = sort(fit.t);
r0 = limbModel(ts,fit.p0);
rf = limbModel(ts,fit.p);

subplot(211)
plot(fit.t,fit.r,'o',ts,r0,'-',ts,rf,'-');
xlabel('\theta [deg]');
if length(fit.p0) == 3,
  title(sprintf('fit: X_c=(%.1f,%.1f) R=%.1f', fit.p))
else
  title(sprintf('fit: X_c=(%.1f,%.1f) a=%.1f b=%.1f tilt=%.0f', fit.p))
end

[legh,objh,oh,om] = legend('data','initial','fitted','location','best');
set(objh(1),'fontsize',9);

% generate model regularly sampled
t = linspace(-180,180,100);
[Sx0,Sy0] = createLimbData(t,length(t),p);

X0 = r0.*cosd(ts);
Y0 = r0.*sind(ts);

Xf = rf.*cosd(ts);
Yf = rf.*sind(ts);

subplot(212)
plot(Sx0,Sy0,'-',Sx,Sy,'o',X0,Y0,'--o',Xf,Yf,'-o');
if length(fit.p0) == 3,
  title(sprintf('model: X_c=(%.1f,%.1f) R=%.1f', p))
else
  title(sprintf('model: X_c=(%.1f,%.1f) a=%.1f b=%.1f tilt=%.0f', p))
end
[legh,objh,oh,om] = legend('model','data+noise','initial','fitted',...
                           'location','best');
set(objh(1),'fontsize',9);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [Sx0,Sy0,varargout] = createLimbData(t,ns,p,pc)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


if length(p) == 3,

  Xc = p(1);
	Yc = p(2);
	r  = p(3);

  % exact circle
	Sx0 = Xc + r*cosd(t);
	Sy0 = Yc + r*sind(t);

  if exist('pc','var') & ~isempty(pc)
    % noisy circle
    Sx = Xc + r*(1+pc*(0.5-rand(1,ns))).*cosd(t);
    Sy = Yc + r*(1+pc*(0.5-rand(1,ns))).*sind(t);
		varargout{1} = Sx;
		varargout{2} = Sy;
	end

elseif length(p) == 5,

  Xc = p(1);
	Yc = p(2);
	a  = p(3);
	b  = p(4);
	t0 = p(5);

	% exact ellipse
	Sx0 = Xc + a*cosd(t)*cosd(t0) - b*sind(t)*sind(t0);
	Sy0 = Yc + a*cosd(t)*sind(t0) + b*sind(t)*cosd(t0);

  if exist('pc','var') & ~isempty(pc)
	  % noisy ellipse
	  Sx = Xc + a*(1+pc*(0.5-rand(1,ns))).*cosd(t)*cosd(t0)-...
	            b*(1+pc*(0.5-rand(1,ns))).*sind(t)*sind(t0);
	  Sy = Yc + a*(1+pc*(0.5-rand(1,ns))).*cosd(t)*sind(t0)+...
	            b*(1+pc*(0.5-rand(1,ns))).*sind(t)*cosd(t0);
		varargout{1} = Sx;
		varargout{2} = Sy;
	end

end

if exist('pc','var') & ~isempty(pc)
  fprintf(1,'%d seeds created with pc %f\n', ns, pc);
  plot(Sx0,Sy0,'x', Sx,Sy,'o')
end
