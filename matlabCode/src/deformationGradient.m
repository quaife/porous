function [time,F11,F12,F21,F22] = deformationGradient(...
    X0,options,prams,fileName)

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

u_x(:,2:end-1) = (u(:,3:end) - u(:,1:end-2))/2/dx;
v_x(:,2:end-1) = (v(:,3:end) - v(:,1:end-2))/2/dx;
u_y(2:end-1,:) = (u(3:end,:) - u(1:end-2,:))/2/dy;
v_y(2:end-1,:) = (v(3:end,:) - v(1:end-2,:))/2/dy;

time = [];
F11 = [];
F12 = [];
F21 = [];
F22 = [];

