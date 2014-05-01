function [Ninner,Nouter,nv,Xinner,Xouter,sigmaInner,sigmaOuter] = ...
    loadFile(fileName)


fid = fopen(fileName,'r');
val = fread(fid,'double');
fclose(fid);

Ninner = val(1);
Nouter = val(2);
nv = val(3);
val = val(4:end);

Xinner = zeros(2*Ninner,nv);
Xouter = zeros(2*Nouter,1);
sigmaInner = zeros(2*Ninner,nv);
sigmaOuter = zeros(2*Nouter,1);

istart = 1;
% start of a pointer to where everything is stored in val

for k = 1:nv
  iend = istart + Ninner - 1;
  Xinner(1:Ninner,k) = val(istart:iend);
  istart = iend + 1;
end
% x-coordinates of the inner boundaries

for k = 1:nv
  iend = istart + Ninner - 1;
  Xinner(Ninner+1:2*Ninner,k) = val(istart:iend);
  istart = iend + 1;
end
% y-coordinates of the inner boundaries

iend = istart + Nouter - 1;
Xouter(1:Nouter,1) = val(istart:iend);
istart = iend + 1;
% x-coordinates of the outer boundary

iend = istart + Nouter - 1;
Xouter(Nouter+1:2*Nouter,1) = val(istart:iend);
istart = iend + 1;
% y-coordinates of the outer boundary


for k = 1:nv
  iend = istart + Ninner - 1;
  sigmaInner(1:Ninner,k) = val(istart:iend);
  istart = iend + 1;
end
% x-coordinates of the inner densities

for k = 1:nv
  iend = istart + Ninner - 1;
  sigmaInner(Ninner+1:2*Ninner,k) = val(istart:iend);
  istart = iend + 1;
end
% y-coordinates of the inner densities 

iend = istart + Nouter - 1;
sigmaOuter(1:Nouter,1) = val(istart:iend);
istart = iend + 1;
% x-coordinates of the outer density

iend = istart + Nouter - 1;
sigmaOuter(Nouter+1:2*Nouter,1) = val(istart:iend);
istart = iend + 1;
% y-coordinates of the outer desnity
