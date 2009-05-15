function [VD, params] = mergeVD(VD, params)
% function [VD, params] = mergeVD(VD, params)

%
% $Id: mergeVD.m,v 1.3 2009/05/15 15:01:36 patrick Exp $
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

% do not attempt to merge if mergePctile < 0
if params.mergePctile<0,
  return;
else
  mergePctile = params.mergePctile;
end

% Similarity parameters |\mu_i-\mu_j|< dmu \mu_i
dmu = params.dmu;
% non homogeneous neighbours vertices length to circumference max ratio
thresHoldLength = params.thresHoldLength;

iMerge = 1;
stopMerge = false;
while ~stopMerge,

  % compute homogeneity fuunction and dynamic threshold
  [WD,SD,WHC,SHC,HCThreshold] = computeHCThreshold(VD, params, mergePctile);

  if 1,
	  [Wmu, VD.Smu] = getVDOp(VD, params.W, @(x) median(x));
  else
	  [Wmu, VD.Smu] = getVDOp(VD, params.W, @(x) mean(x));
  end

  if 0, % diagnostic plot
    [vx,vy] = voronoi(VD.Sx(VD.Sk), VD.Sy(VD.Sk));
  end

  % build index table to map values in VD.Nk to their indices 
  IST = zeros(length(VD.Nk),1);
  i=1;
  for in = 1:length(VD.Nk),
    if ~isempty(find(VD.Sk == in)),
	    IST(in) = i;
	    i = i+1;
	  end
  end

  Sk = [];
  for isk = IST(VD.Sk(find(SHC < HCThreshold)))', % For homogeneous VR
    % isk contains indice in array of size number of seeds
	  %  sk contains seed indice in VD.Nk
    sk = VD.Sk(isk);
    if 0,
	    fprintf(1,'s=%3d HC=%5.2f (%d) mu=%5.2f HC Threshold=%5.2f\n', ...
		          [sk, SHC(isk),(SHC(isk) < HCThreshold),VD.Smu(isk),HCThreshold]);
	  end
	  % Homogeneous neighbour VR
	  ihc = (SHC(IST(VD.Nk{sk}))' < HCThreshold);
	  if 0,
	    fprintf(1,'n=%3d HC=%5.2f (%d) mu=%5.2f\n', ...
	            [VD.Nk{sk}, SHC(IST(VD.Nk{sk})), ihc', VD.Smu(IST(VD.Nk{sk}))]');
	  end
	  % 
    err  =  abs(VD.Smu(isk) - VD.Smu(IST(VD.Nk{sk}(ihc)))'); 
	  if 0,
	    if ~isempty(err), 
	      fprintf(1,'err=');
	      fprintf(1,' %5.2f', err); 
        fprintf(1,'\n');
	    end
	  end
	  if 0, % diagnostic plot
	    subplot(111)
	    imagesc(Wmu), 
	    set(gca,'clim',params.Wlim);
	    colorbar
	    axis xy,
	    hold on
	    plot(VD.Sx(sk), VD.Sy(sk), 'xk', 'MarkerSize', 5);
	    plot(VD.Sx(VD.Nk{sk}(ihc)), VD.Sy(VD.Nk{sk}(ihc)), 'ok', 'MarkerSize', 5);
	    plot(vx,vy,'-k','LineWidth',0.5)
	    hold off
	    %pause
	  end
	  if all(err < dmu*abs(VD.Smu(isk))), 
	    % all homogeneous neighbours are less than dmu\% different
      % get vertices list of VR(sk)
		  % unbounded VR not handled properly
      [V,I] = getVRvertices(VD, sk);
		  edgesLength = sqrt(sum((V([1:end],:)-V([2:end 1],:)).^2,2));
			if 0,
		  fprintf(1,'edgesLength='); 
			fprintf(1,' %5.2f', edgesLength);
			fprintf(1,'\n');
			end
	  	totalLength = sum(edgesLength);
		  nonHCLength = sum(edgesLength(I(ihc)));
		  r = nonHCLength/totalLength;
			if 0
		  fprintf(1,'totalLength=%.2f nonHCLength=%.2f ratio=%f\n',...
			        totalLength, nonHCLength, r);
			end
		  if r < thresHoldLength,
			  if 0
		    fprintf(1,'s=%d to be removed\n', sk);
				end
		    Sk = [Sk; sk]; 
		  end
		  %pause
	  end
  end
  %save([params.oDir 'removedSeeds'], 'Sk');

	% save homogeneity function and dynamic threshold
	VD.mergeSHC{iMerge} = SHC;
	VD.mergeHCThreshold(iMerge) = HCThreshold;
  if ~isempty(Sk),
    fprintf(1,'Removing %d seeds to Voronoi Diagram\n', length(Sk));
    %pause
    for k = Sk(:)',
      VD  = removeSeedFromVD(VD, k);
			fprintf(1,'Voronoi Diagram computed\n')
	    if 0
        drawVD(VD);
	    end
    end
    params = plotCurrentVD(VD, params, iMerge);
	  iMerge = iMerge+1;
  else
    stopMerge = true;
  end
end


function params = plotCurrentVD(VD, params, iMerge)

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
plot(vx,vy,'-k','LineWidth',1)
hold off

title(sprintf('card(S) = %d  (iteration %d)', length(VD.Sk), iMerge))

drawnow 

if params.mergeExport,
  exportfig(gcf,[params.oDir 'merge' num2str(iMerge) '.eps'],'color','cmyk');
end

if params.movDiag,
  params.mov = addframe(params.mov, getframe(gcf,[0 0 params.movPos(3:4)]));
end

