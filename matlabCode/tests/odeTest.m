x0 = -2*ones(40,1);
y0 = linspace(-.2,.2,numel(x0))';
%%y0 = linspace(-1/2,1/2,numel(x0))';
%x0 = -2; y0 = 0.3;
%x0 = -3.581193390559622e-01;
%y0 = 9.336758216667427e-01;
X0 = [x0;y0];


T = 20;
nT = 201;
options.RelTol = 1e-3;
options.AbsTol = 1e-6;

[t,X] = ode45(@exactODEfun,linspace(0,T,nT),X0,options);
xExact = X(:,1:end/2); yExact = X(:,end/2+1:end);
% using exact velocity field

[t,X] = ode45(@cartInterpODEfun,linspace(0,T,nT),X0,options);
xCart = X(:,1:end/2); yCart = X(:,end/2+1:end);
% using an interpolation of a cartesian grid

[t,X] = ode45(@radInterpODEfun,linspace(0,T,nT),X0,options);
xRad = X(:,1:end/2); yRad = X(:,end/2+1:end);
% using an interpolation of a polar grid

[r,theta] = meshgrid(1:.1:3,0:2*pi/20:2*pi);
r = r(:); theta = theta(:);
xgrid = r.*cos(theta); ygrid = r.*sin(theta);
vel = exactODEfun(0,[xgrid;ygrid]);
u = vel(1:end/2); v = vel(end/2+1:end);

figure(1);
clf; hold on
quiver(xgrid,ygrid,u,v);
theta = linspace(0,2*pi,100);
plot(cos(theta),sin(theta),'k')

for k = 1:numel(x0)
  plot(xExact(:,k),yExact(:,k),'r');
end
axis equal
axis(1.5*[-1 1 -1 1])
%axis(ax)



