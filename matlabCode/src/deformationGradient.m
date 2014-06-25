%function [time,F11,F12,F21,F22] = deformationGradient(...
%    X0,options,prams,fileName)

om = monitor(options,prams);

[Ninner,Nouter,nv,Xinner,Xouter,sigmaInner,sigmaOuter] = ...
    om.loadGeometry(fileName);

fileName1 = [fileName(1:end-8) 'EulerVelocities.bin'];
[ny,nx,eulerX,eulerY,u,v] = om.loadEulerVelocities(fileName1);

u_x = zeros(ny,nx);
u_y = zeros(ny,nx);
v_x = zeros(ny,nx);
v_y = zeros(ny,nx);

dx = eulerX(1,2) - eulerX(1,1);
dy = eulerY(2,1) - eulerY(1,1);

u_x(:,1:end-1) = (u(:,2:end) - u(:,1:end-1))/dx;
u_y(1:end-1,:) = (u(2:end,:) - u(1:end-1,:))/dy;
v_x(:,1:end-1) = (v(:,2:end) - v(:,1:end-1))/dx;
v_y(1:end-1,:) = (v(2:end,:) - v(1:end-1,:))/dy;
% first-order finite difference for every point except the final
% column/row

u_x(:,end) = (u(:,end) - u(:,end-1))/dx;
v_x(:,end) = (u(:,end) - u(:,end-1))/dx;
u_y(end,:) = (u(end,:) - u(end-1,:))/dy;
v_y(end,:) = (v(end,:) - v(end-1,:))/dy;
% deal with the boundaries by using a finite difference going in the
% other direction

indy = 4835; indx = 234;
x = eulerX(indy,indx);
y = eulerY(indy,indx);
% Some random point that is well separated from the pores

time = linspace(0,prams.T,prams.ntime);
dt = time(2) - time(1);
F11 = zeros(size(time));
F21 = zeros(size(time));
F12 = zeros(size(time));
F22 = zeros(size(time));

F11(1) = 1;
F22(1) = 1;
% initial condition is identity matrix

op = poten(0,false);
Jacobian = [u_x(indy,indx);u_y(indy,indx);...
            v_x(indy,indx);v_y(indy,indx)];
odeFun = @(t,z) op.deformationGradientRHS(t,z,Jacobian);
F0 = [1;0;0;1];
[time,F] = ode45(odeFun,time,F0);

%for k = 1:numel(time)-1
%  F11(k+1) = F11(k) + dt*(...
%      u_x(indy,indx)*F11(k) + u_y(indy,indx)*F21(k));
%  F21(k+1) = F21(k) + dt*(...
%      u_x(indy,indx)*F12(k) + u_y(indy,indx)*F22(k));
%  F12(k+1) = F12(k) + dt*(...
%      v_x(indy,indx)*F11(k) + v_y(indy,indx)*F21(k));
%  F22(k+1) = F22(k) + dt*(...
%      v_x(indy,indx)*F12(k) + v_y(indy,indx)*F22(k));
%end
%




%theta = (0:31)*2*pi/32;
%alpha = 1e-1;
%xbody = x + alpha*cos(theta);
%ybody = y + alpha*sin(theta);
%clf
%plot(xbody,ybody,'r')
%axis equal
%pause
%
%for k = 1:numel(time)-1
%  z = [xbody;ybody];
%  F = [F11(k) F12(k); F21(k) F22(k)]; 
%  Deltaz = F*z;
%  z = z + dt*Deltaz;
%  xbody = z(1,:);
%  ybody = z(2,:);
%  plot(xbody,ybody,'r')
%  axis equal
%  pause
%
%end




