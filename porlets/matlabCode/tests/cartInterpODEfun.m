function fX = cartInterpODEfun(t,X)

x = X(1:end/2); y = X(end/2+1:end);

n = 128;
[eX,eY] = meshgrid(linspace(-2,2,n),linspace(-2,2,n));
u1 = -1/4*log(eX.^2 + eY.^2);
u2 = (2*eX.^4 + 2*eX.^2.*eY.^2 - eX.^2 + eY.^2)/4./(eX.^2+eY.^2).^2;
v1 = 0;
v2 = eX.*eY.*(eX.^2+eY.^2-1)/2./(eX.^2+eY.^2).^2;
uinf = 0.25;
vinf = 0.0;

eU = uinf - (u1+u2);
eV = vinf - (v1+v2);

s = find(eX.^2 + eY.^2 <=1);
eU(s) = 0;
eV(s) = 0;

u = interp2(eX,eY,eU,x,y,'spline');
v = interp2(eX,eY,eV,x,y,'spline');


fX = [u;v];





