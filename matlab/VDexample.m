function VDexample
% function VDexample

% init seed of Mersenne-Twister RNG
if exist('RandStream','file') == 2,
  RandStream.setDefaultStream(RandStream('mt19937ar','seed',30)));
else
  rand('twister',30);
end

x = round(100*rand(8,1))+1;
y = round(100*rand(8,1))+1;

xm = min(x); xM = max(x);
ym = min(y); yM = max(y);


[vx,vy] = voronoi(x, y);

VD = computeVD(100,100, [x, y]);

hold on
plot(x,y,'xk','MarkerSize',5)
plot(vx,vy,'-k','LineWidth',1);
set(gca,'xlim',[-1 2],'ylim',[-1 2])
for i=1:length(x)
  text(x(i), y(i), num2str(i), 'verticalalignment', 'bottom');
  for j=setdiff(VD.Nk{i}', 1:i),
    plot([x(i), x(j)], [y(i), y(j)], ':k','LineWidth',.5)
  end
end
hold off
axis off
set(gca,'xlim',100*[-.5 1.5],'ylim',100*[-.5 1.5]);
axis equal

%orient tall
%exportfig(gcf,'../agu2008/fig1.eps','color','cmyk');
%opts = struct('color','cmyk','bounds','tight');
%opts = struct('color','cmyk','bounds','tight','linestylemap','bw');
%exportfig(gcf,'../agu2008/fig1.eps',opts);
%savefig('fig1.eps',gcf,'eps','-cmyk','-crop');
