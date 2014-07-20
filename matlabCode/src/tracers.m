function [time,xtra,ytra,F11,F12,F21,F22] = tracers(...
    X0,options,prams,fileName)

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
nparts = options.nparts;
% size of domain for Eulerian grid and the number of sections the
% Eulerian grid is divided into.  If it is too large, memory becomes a
% problem

if options.computeEuler
  [eulerX,eulerY] = ...
      meshgrid(linspace(xmin,xmax,nx),linspace(ymin,ymax,ny));
  % build the Eulerian grid that is used to interpolate in the time
  % integrator

  cutoff = ceil(numel(eulerX)/nparts);
  % maximum number of points to compute at once
  eX = eulerX(:); eY = eulerY(:);

  tic
  vel = zeros(2*numel(eX),1);
  for k = 1:nparts
    istart = (k-1)*cutoff + 1;
    iend = min(istart + cutoff - 1,numel(eX));
    velPart = op.layerEval(0,[eX(istart:iend);eY(istart:iend)],...
        options.ymThresh,options.ypThresh,...
        innerGeom,outerGeom,sigmaInner,sigmaOuter);
    % compute velocity at points [eX;eY]
    vel(istart:iend) = velPart(1:end/2);
    vel((istart:iend)+numel(eX)) = velPart(end/2+1:end);
  end
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
  % load velocities at Eulerian grid from a precomputed computation
  if (nx ~= options.nx || ny ~= options.ny)
    message = 'Saved Euler grid does not match input parameters';
    om.writeMessage(message);
    message = 'Just an FYI';
    om.writeMessage(message);
  end
end


[u_x,u_y,v_x,v_y] = computeDerivs(eulerX,eulerY,u,v);
% find gradient of velocity field using finite differences.  Want to use
% one-sided derivatives on the boundary and when the stencil is cut by a
% pore


%figure(1); clf;
%contourf(eulerX(1:20:end,1:20:end),eulerY(1:20:end,1:20:end),log10(abs(u_x(1:20:end,1:20:end) + v_y(1:20:end,1:20:end))))
%colorbar
%
%[dx dy]
%max(max(abs(u_x + v_y)))
%max(max(abs(u_x(2:end-1,2:end-1) + v_y(2:end-1,2:end-1))))
%pause

odeFun = @(t,z) op.interpolateLayerPot(t,z,eulerX,eulerY,...
    u,v,u_x,u_y,v_x,v_y,prams.T,options.ymThresh);
% function handle that evalutes the right-hand side.  Handles the
% position and deformation gradient all at once.

tic
opts.RelTol = prams.rtol;
opts.AbsTol = prams.atol;
% load options for ode45

ntra = numel(X0)/2;
xtra = zeros(prams.ntime,ntra);
ytra = zeros(prams.ntime,ntra);
F11 = zeros(prams.ntime,ntra);
F12 = zeros(prams.ntime,ntra);
F21 = zeros(prams.ntime,ntra);
F22 = zeros(prams.ntime,ntra);
% allocate memory for positions and deformation gradient

for k = 1:numel(X0)/2
  message = ['\ntracers ' num2str(2*k/numel(X0)*100,'%04.1f\n') ' %% completed\n'];
  fprintf(message);
  x0 = X0(k);
  y0 = X0(k+numel(X0)/2);
  % initial condition of the kth tracer

  [time,z] = ode45(odeFun,linspace(0,prams.T,prams.ntime),...
      [x0 y0 1 0 0 1],opts);
  % Solve for position and deformation gradient with ode45

  xtra(:,k) = z(:,1);
  ytra(:,k) = z(:,2);
  F11(:,k) = z(:,3);
  F12(:,k) = z(:,4);
  F21(:,k) = z(:,5);
  F22(:,k) = z(:,6);

%  figure(1); clf;hold on
%  plot(Xouter(1:end/2),Xouter(end/2+1:end),'k')
%  fill(Xinner(1:end/2,:),Xinner(end/2+1:end,:),'k')
%  plot(xtra(:,1:k),ytra(:,1:k),'r-')
%  axis equal
%  ax = [min(xtra(:,k))-0.2 max(xtra(:,k))+0.2 ...
%      min(ytra(:,k))-0.2 max(ytra(:,k))+0.2];
%  axis(ax)
%
%  figure(2); clf; hold on
%  plot(F11(:,k))
%  plot(F12(:,k))
%  plot(F21(:,k))
%  plot(F22(:,k))

  if mod(k,1000) == 1
    om.writeTracerPositions(time,xtra(:,1:k),ytra(:,1:k));
    om.writeDeformationGradient(time,F11(:,1:k),F12(:,1:k),...
        F21(:,1:k),F22(:,1:k));
  end
  % save every 100th iteration
end
om.writeMessage(' ');

om.writeStars
message = '****       Tracer locations found        ****';
om.writeMessage(message);
message = ['**** Required time was ' num2str(toc,'%4.2e') ...
    ' seconds  ****'];
om.writeMessage(message);
om.writeStars
om.writeMessage(' ');

if options.usePlot
  om.plotData;
  om.runMovie;
end

om.writeTracerPositions(time,xtra,ytra);
om.writeDeformationGradient(time,F11,F12,F21,F22);
% one final save at the very end



