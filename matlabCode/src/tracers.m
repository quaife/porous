function [t,xtra,ytra] = tracers(X0,options,prams,fileName)

om = monitor(options,prams);

[Ninner,Nouter,nv,Xinner,Xouter,sigmaInner,sigmaOuter] = ...
    loadFile(fileName);
% load all the information about the geometry and density function
if ((prams.Ninner == Ninner) + ...
      (prams.Nouter == Nouter) + (prams.nv == nv))~=3
  message = 'Saved data does not match input parameters';
  om.writeMessage(message);
  message = 'I am stopping';
  om.writeMessage(message);
  return
end

innerGeom = capsules(Xinner,'inner');
outerGeom = capsules(Xouter,'inner');
% build objects for the inner and outer boundaries

fmm = options.fmm;
op = poten(Ninner,fmm);

xmin = options.xmin; xmax = options.xmax; nx = options.nx;
ymin = options.ymin; ymax = options.ymax; ny = options.ny;
[eulerX,eulerY] = ...
    meshgrid(linspace(xmin,xmax,nx),linspace(ymin,ymax,ny));
% build the Eulerian grid that is used to interpolate in the time
% integrator

tic
vel = op.layerEval(0,[eulerX(:);eulerY(:)],...
    options.ymThresh,options.ypThresh,...
    innerGeom,outerGeom,sigmaInner,sigmaOuter);
% evalute velocity on an Eulerian grid
om.writeStars
message = '****   Velocity found on Eulerian Grid   ****';
om.writeMessage(message);
message = ['**** Required time was ' num2str(toc,'%4.2e') ...
    ' seconds  ****'];
om.writeMessage(message);
om.writeStars
om.writeMessage(' ');

u = reshape(vel(1:end/2),size(eulerX));
v = reshape(vel(end/2+1:end),size(eulerY));
% put velocity field in format that works well for interp2
%figure(1); clf; hold on
%plot(Xinner(1:end/2,:),Xinner(end/2+1:end,:),'k')
%quiver(eulerX,eulerY,u,v,'r')
%pause

odeFun = @(t,z) op.interpolateLayerPot(t,z,eulerX,eulerY,u,v);
% function handle that evalutes the right-hand side 
tic
opts.RelTol = 1e-5;
opts.AbsTol = 1e-8;
[t,Xtra] = ode45(odeFun,linspace(0,20,200),X0);

om.writeStars
message = '****       Tracer locations found        ****';
om.writeMessage(message);
message = ['**** Required time was ' num2str(toc,'%4.2e') ...
    ' seconds  ****'];
om.writeMessage(message);
om.writeStars
om.writeMessage(' ');

xtra = Xtra(:,1:end/2);
ytra = Xtra(:,end/2+1:end);

