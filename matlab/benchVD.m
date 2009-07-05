function benchVD
% function benchVD

%
% $Id: benchVD.m,v 1.3 2009/07/05 19:44:03 patrick Exp $
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

%ns = fix(logspace(log10(10),log10(3000),30));
ns = fix(linspace(10,3000,30));
%ns = fix(linspace(20,300,30));
nr = 256;
nc = 256;
initSeeds = @randomSeeds;

% init seed of Mersenne-Twister RNG
rand('twister',10);

if exist('initSeeds') & isa(initSeeds, 'function_handle'),
  [initSeeds, msg] = fcnchk(initSeeds);
  S = initSeeds(nr, nc, ns(end));
else
  error('initSeeds not defined or not a Function Handle');
end

nd = [ns(1), diff(ns)];

tVDa = zeros(size(ns));
tVDf = zeros(size(ns));

% initial seeds
s = S([1:ns(1)],:);

% incremental add
tStart = tic;
VDa = computeVD(nr, nc, s);
tVDa(1) = toc(tStart);
fprintf(1,'add    %4d seeds (%3d:%3d) %8.1f s\n', size(s,1), 1, ns(1), tVDa(1));

% full 
tStart = tic;
VDf = computeVDFast(nr, nc, s);
tVDf(1) = toc(tStart);
fprintf(1,'full   %4d seeds (%3d:%3d) %8.1f s\n', size(s,1), 1, ns(1), tVDf(1));

for i=2:length(ns)

  % incremental add
  s = S([ns(i-1)+1:ns(i)],:);
  tStart = tic;
  for k = 1:length(s),
    VDa = addSeedToVD(VDa, s(k,:));
	end
  tVDa(i) = toc(tStart);
  fprintf(1,'add    %4d seeds (%3d:%3d) %8.1f s\n', ...
          size(s,1),ns(i-1)+1,ns(i),tVDa(i));

  % full
  s = S([1:ns(i)],:);
  tStart = tic;
  VDf = computeVDFast(nr, nc, s);
  tVDf(i) = toc(tStart);
  fprintf(1,'full   %4d seeds (%3d:%3d) %8.1f s\n', size(s,1), 1, ns(i), tVDf(i));

  plot(ns(1:i), tVDa(1:i)./nd(1:i), ns(1:i), tVDf(1:i));
	drawnow

end

% incremental remove
nsr = ns(end:-1:1);
ndr = -diff(nsr);

tVDr = zeros(size(nsr));

VDr = VDa;
for i=1:length(nsr)-1,
  
  s = VDr.Sk(nsr(i):-1:nsr(i+1)+1);
  tStart = tic;
  for k = 1:length(s),
    VDr = removeSeedFromVD(VDr, s(k,:));
	end
  tVDr(i) = toc(tStart);

  fprintf(1,'remove %4d seeds (%3d:%3d) %8.1f s\n', ...
	        size(s,1), nsr(i), nsr(i+1)+1, tVDr(i));

  plot(nsr(1:i), tVDr(1:i)./ndr(1:i));
	drawnow
end
nsr(end) = [];
tVDr(end) = [];

[ptVDa] = polyfit(ns, tVDa./nd, 1);
[ptVDf] = polyfit(ns, tVDf, 2);
[ptVDr] = polyfit(nsr, tVDr./ndr, 1);

subplot(211),
plot(ns, [tVDf; polyval(ptVDf,ns)], '-o');
xlabel('number of seeds')
ylabel('time [s]')
title('Full VOISE')

subplot(212),
plot(ns, [tVDa./nd; polyval(ptVDa,ns)], '-o', ...
		 nsr, [tVDr./ndr; polyval(ptVDr,nsr)], '-o');
legend('Add','Add fit','Remove','Remove fit','location','northwest')
xlabel('number of seeds')
ylabel('time [s]')
title('Incremental VOISE')


save VOISEtiming ns nd tVDa tVDf ptVDa ptVDf nsr ndr tVDr ptVDr



