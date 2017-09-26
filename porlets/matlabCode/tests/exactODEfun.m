function fX = exactODEfun(t,X)

x = X(1:end/2); y = X(end/2+1:end);
u1 = -1/4*log(x.^2 + y.^2);
u2 = (2*x.^4 + 2*x.^2.*y.^2 - x.^2 + y.^2)/4./(x.^2+y.^2).^2;
v1 = 0;
v2 = x.*y.*(x.^2+y.^2-1)/2./(x.^2+y.^2).^2;
uinf = 0.25;
vinf = 0.0;

fX = [uinf - (u1+u2);vinf - (v1+v2)];


