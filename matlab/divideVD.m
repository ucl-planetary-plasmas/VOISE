function [VD, params]  = divideVD(VD, params)
% function [VD,params] = divideVD(VD, params)

%
% $Id: divideVD.m,v 1.2 2009/03/18 15:45:56 patrick Exp $
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

% do not attempt to divide if dividePctile < 0
if params.dividePctile<0,
  return;
else
  dividePctile = params.dividePctile;
end

iDiv = 1;
stopDiv = false;
while ~stopDiv,

  % compute homogeneity fuunction and dynamic threshold
  [WD,SD,WHC,SHC,HCThreshold] = computeHCThreshold(VD, params, dividePctile);

  if 0, % diagnostic plot (histogram)
    subplot(211),
    edges = linspace(min(SHC),max(SHC),10);
    n = histc(SHC,edges); bar(edges, cumsum(n)/sum(n)*100, 'histc')
    %pause
	end

  S = ones(0,2);
  for sk = VD.Sk(find(SHC >= HCThreshold))', 
	  % check only the nonhomogeneous Voronoi regions
    s = addSeedsToVR(VD, sk, params);
    S = [[S(:,1); s(:,1)], [S(:,2); s(:,2)]];
  end
	
	if 0, % diagnostic plot (seed to add)
    subplot(212),
    imagesc(WHC);
    axis xy, colorbar
    hold on
    for i=1:size(S,1),
      plot(S(i,1), S(i,2), 'ok', 'MarkerSize', 5);
    end
    hold off
    %pause
	end

	% save homogeneity function and dynamic threshold
	VD.divSHC{iDiv} = SHC;
	VD.divHCThreshold(iDiv) = HCThreshold;
  if ~isempty(S),
    fprintf(1,'Adding %d seeds to Voronoi Diagram\n', size(S,1))
    for k = 1:size(S,1),
		  if isempty(find(S(k,1)==VD.Sx & S(k,2)==VD.Sy)),
        VD = addSeedToVD(VD, S(k,:));
				if 0, % diagnostic plot
          drawVD(VD);
				end
			end
    end
	  params = plotCurrentVD(VD, params, iDiv);
	  iDiv = iDiv+1;
  else
    stopDiv = true;
  end
  fprintf(1,'Voronoi Diagram computed\n')
	if 0, % diagnostic plot
    drawVD(VD);
    %pause
	end
end

function params = plotCurrentVD(VD, params, iDiv)

VDW = getVDOp(VD, params.W, @(x) median(x));

clf
subplot(111),
imagesc(VDW),
axis xy,
axis equal
axis off
set(gca,'clim',params.Wlim);
%colorbar
set(gca,'xlim',[VD.xm VD.xM], 'ylim', [VD.ym VD.yM]);

hold on
[vx,vy]=voronoi(VD.Sx(VD.Sk), VD.Sy(VD.Sk));
plot(vx,vy,'-k','LineWidth',0.5)
hold off

title(sprintf('card(S) = %d  (iteration %d)', length(VD.Sk), iDiv))

drawnow
 
if params.divideExport,
  exportfig(gcf,[params.oDir 'div' num2str(iDiv) '.eps'], 'color','cmyk');
end

if params.movDiag,
  params.mov = addframe(params.mov, getframe(gcf,[0 0 params.movPos(3:4)]));
end

