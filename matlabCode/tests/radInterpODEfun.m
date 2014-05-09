function fX = cartInterpODEfun(t,X)

x = X(1:end/2); y = X(end/2+1:end);

n = 4;
[r,theta] = meshgrid(linspace(1,sqrt(8),n),linspace(0,2*pi*(1+1/n),n));
% add some overlap for periodicity
eX = r.*cos(theta); eY = r.*sin(theta);

u1 = -1/4*log(eX.^2 + eY.^2);
u2 = (2*eX.^4 + 2*eX.^2.*eY.^2 - eX.^2 + eY.^2)/4./(eX.^2+eY.^2).^2;
v1 = 0;
v2 = eX.*eY.*(eX.^2+eY.^2-1)/2./(eX.^2+eY.^2).^2;
uinf = 0.25;
vinf = 0.0;

eU = uinf - (u1+u2);
eV = vinf - (v1+v2);

z = x + 1i*y;
r0 = abs(z);
theta0 = angle(z);
theta0(theta0 < 0) = 2*pi+theta0(theta0<0);

u = interp2(r,theta,eU,r0,theta0);
v = interp2(r,theta,eV,r0,theta0);

fX = [u;v];





