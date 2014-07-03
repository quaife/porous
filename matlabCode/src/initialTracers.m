function [xtar,ytar] = initialTracers(radii,centers,ntra)
% find a set of intial tracer locations that are in the middle of the
% porous region
% keep xtar between 0.5 and 4.3
% keep ytar between 5 and 30

xmin = 0.5; xmax = 4.3; dx = xmax - xmin;
ymin = 5; ymax = 30; dy = ymax - ymin;
xtar = [];
ytar = [];


while numel(xtar) < ntra
  xtarPot = dx*rand(ceil(1.1*(ntra-numel(xtar))),1) + xmin;
  ytarPot = dy*rand(ceil(1.1*(ntra-numel(ytar))),1) + xmax;

  ncount = 0;
  nv = numel(radii);
  for k = 1:numel(xtarPot)
    if(~any((xtarPot(k) - centers(1:nv,1)).^2 + ...
          (ytarPot(k) - centers(1:nv,2)).^2 < radii(1:nv).^2))
      xtar = [xtar xtarPot(k)]; 
      ytar = [ytar ytarPot(k)]; 
    end
  end
end

xtar = xtar(1:ntra);
ytar = ytar(1:ntra);


