function benchVD
% function benchVD

%
% $Id: benchVD.m,v 1.2 2009/07/04 15:04:16 patrick Exp $
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

%is = [10 20 40 80 160 320 640 1280 2560];
is = fix(logspace(log10(10),log10(3000),30));
nr = 256;
nc = 256;
ns = is(end);
initSeeds = @randomSeeds;

% init seed of Mersenne-Twister RNG
rand('twister',10);

if exist('initSeeds') & isa(initSeeds, 'function_handle'),
  [initSeeds, msg] = fcnchk(initSeeds);
  S = initSeeds(nr, nc, ns);
else
  error('initSeeds not defined or not a Function Handle');
end



tVDi = zeros(size(is));
tVDf = zeros(size(is));

s = S([1:is(1)],:);

% incremental 
tStart = tic;
VD = computeVD(nr, nc, s);
tVDi(1) = toc(tStart);

% full
tStart = tic;
VDf = computeVDFast(nr, nc, s);
tVDf(1) = toc(tStart);

fprintf(1,'%4d seeds: elapsed time inc/full %8.3f/%8.3f\n', ...
        size(s,1), tVDi(1), tVDf(1));

for i=2:length(is)

  % incremental
  s = S([is(i-1)+1:is(i)],:);
  tStart = tic;
  for k = 1:length(s),
    VD = addSeedToVD(VD, s(k,:));
	end
  tVDi(i) = toc(tStart);

  % full
  s = S([1:is(i)],:);
  tStart = tic;
  VDf = computeVDFast(nr, nc, s);
  tVDf(i) = toc(tStart);

  fprintf(1,'%4d seeds: elapsed time inc/full %8.3f/%8.3f\n', ...
	        size(s,1), tVDi(i), tVDf(i));

  plot(is(1:i), tVDi(1:i), is(1:i), tVDf(1:i));
	drawnow

end


ptVDi = polyfit(is, tVDi, 2);
ptVDf = polyfit(is, tVDf, 2);

plot(is, [tVDi;polyval(ptVDi,is)], '-o', is, [tVDf;polyval(ptVDf,is)], '-o');
xlabel('number of seeds')
ylabel('time [s]')

save VOISEtiming is tVDi tVDf ptVDi ptVDf


