function [time,xtra,ytra] = tracers(X0,options,prams,fileName)

om = monitor(options,prams);

[Ninner,Nouter,nv,Xinner,Xouter,sigmaInner,sigmaOuter] = ...
    om.loadGeometry(fileName);
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
outerGeom = capsules(Xouter,'outer');
% build objects for the inner and outer boundaries

fmm = options.fmm;
op = poten(Ninner,fmm);

xmin = options.xmin; xmax = options.xmax; nx = options.nx;
ymin = options.ymin; ymax = options.ymax; ny = options.ny;

if options.computeEuler
  [eulerX,eulerY] = ...
      meshgrid(linspace(xmin,xmax,nx),linspace(ymin,ymax,ny));
  % build the Eulerian grid that is used to interpolate in the time
  % integrator

%  [r,theta] = ...
%    meshgrid(linspace(1.37607e-1,1.0e0,100),(0:99)*2*pi/100);
%  eulerX = r.*cos(theta) + 3.9010990e0;
%  eulerY = r.*sin(theta) + 2.4065934e1;

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
  om.writeEulerVelocities(eulerX,eulerY,u,v);
  % save the velocity field so we don't have to keep recomputing it
else
  fileName1 = [fileName(1:end-8) 'EulerVelocities.bin'];
  [ny,nx,eulerX,eulerY,u,v] = om.loadEulerVelocities(fileName1);
  if (nx ~= options.nx || ny ~= options.ny)
    message = 'Saved Euler grid does not match input parameters';
    om.writeMessage(message);
    message = 'Just an FYI';
    om.writeMessage(message);
  end
end


odeFun = @(t,z) op.interpolateLayerPot(t,z,eulerX,eulerY,u,v,prams.T);
% function handle that evalutes the right-hand side 
tic
opts.RelTol = prams.rtol;
opts.AbsTol = prams.atol;
[time,Xtra] = ode45(odeFun,linspace(0,prams.T,prams.ntime),X0);
om.writeMessage(' ');

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

if options.usePlot
%  om.plotData(Xinner,Xouter,eulerX,eulerY,u,v,xtra,ytra);
%  om.runMovie(Xinner,Xouter,xtra,ytra);
  om.plotData;
  om.runMovie;
end

om.writeTracerPositions(time,xtra,ytra);
fileName1 = [fileName(1:end-8) 'TracerPositions.bin'];
[ntime,ntra,time,xtra,ytra] = ...
    om.loadTracerPositions(fileName1);


