function [] = tracers(X0,options,prams,radii,centers,fileName)
%function [time,xtra,ytra,F11,F12,F21,F22] = tracers(...
%    X0,options,prams,fileName)

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
op = poten(innerGeom,fmm,false,options.computeEuler);

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

  vel = zeros(2*numel(eX),1);
%  [nx,ny,eulerX,eulerY,u,v] = om.loadEulerVelocities('/scratch/quaife/porousSimulations/results/newGeoms/circles37EulerVelocities.bin');
%  eX = eulerX(:); eY = eulerY(:);
%  vel = [u(:);v(:)];
% can be used to restart a simulation part way through

  tic
  for k = 1:50
    disp([k nparts])
    istart = (k-1)*cutoff + 1;
    iend = min(istart + cutoff - 1,numel(eX));

    velPart = op.layerEval(0,[eX(istart:iend);eY(istart:iend)],...
        options.xmThresh,options.xpThresh,...
        innerGeom,outerGeom,sigmaInner,sigmaOuter);
    % compute velocity at points [eX;eY]
    vel(istart:iend) = velPart(1:end/2);
    vel((istart:iend)+numel(eX)) = velPart(end/2+1:end);

    u = reshape(vel(1:end/2),size(eulerX));
    v = reshape(vel(end/2+1:end),size(eulerY));
   % put velocity field in format that works well for interp2
    om.writeEulerVelocities(eulerX,eulerY,u,v);
    % save the velocity field so we don't have to keep recomputing it
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
  if options.defGradient
    [u_x,u_y,v_x,v_y] = computeDerivs(eulerX,eulerY,u,v);
    % find gradient of velocity field using finite differences.  Want to
    % use one-sided derivatives on the boundary and when the stencil is
    % cut by a pore
    om.writeEulerDerivs(u_x,u_y,v_x,v_y);
    % save the gradient of the velocity field
  else
    u_x = []; u_y = []; v_x = []; v_y = [];
  end
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
  if options.defGradient
    fileName1 = [fileName(1:end-8) 'EulerDerivs.bin'];
    [ny,nx,u_x,u_y,v_x,v_y] = om.loadEulerDerivs(fileName1);
    % load gradient of velocity on Eulerian grid
    if (nx ~= options.nx || ny ~= options.ny)
      message = 'Saved Euler grid does not match input parameters';
      om.writeMessage(message);
      message = 'Just an FYI';
      om.writeMessage(message);
    end
  else
    u_x = []; u_y = []; v_x = []; v_y = [];
  end
end

odeFun = @(t,z) op.interpolateLayerPot(t,z,eulerX,eulerY,...
    u,v,u_x,u_y,v_x,v_y,prams.T,...
    options.xmThresh,options.xpThresh,...
    options.ymThresh,options.ypThresh,...
    options.defGradient);
% function handle that evalutes the right-hand side.  Handles the
% position and deformation gradient all at once.

xWin = eulerX(eulerX > options.xmThresh + 0.1 & ...
    eulerX < options.xmThresh + 2 & ...
    eulerY > options.ymThresh + 0.1 & ...
    eulerY < options.ypThresh - 0.1);
yWin = eulerY(eulerX > options.xmThresh + 0.1 & ...
    eulerX < options.xmThresh + 2 & ...
    eulerY > options.ymThresh + 0.1 & ...
    eulerY < options.ypThresh - 0.1);
uWin = u(eulerX > options.xmThresh + 0.1 & ...
    eulerX < options.xmThresh + 2 & ...
    eulerY > options.ymThresh + 0.1 & ...
    eulerY < options.ypThresh - 0.1);
vWin = v(eulerX > options.xmThresh + 0.1 & ...
    eulerX < options.xmThresh + 2 & ...
    eulerY > options.ymThresh + 0.1 & ...
    eulerY < options.ypThresh - 0.1);
xWin = xWin(:); yWin = yWin(:); uWin = uWin(:); vWin = vWin(:);

tic
opts.RelTol = prams.rtol;
opts.AbsTol = prams.atol;
% load options for ode45

ntra = numel(X0)/2;
xtraLinear = zeros(prams.ntime,ntra);
ytraLinear = zeros(prams.ntime,ntra);
F11Linear = zeros(prams.ntime,ntra);
F12Linear = zeros(prams.ntime,ntra);
F21Linear = zeros(prams.ntime,ntra);
F22Linear = zeros(prams.ntime,ntra);
xtraLog = zeros(prams.ntime,ntra);
ytraLog = zeros(prams.ntime,ntra);
F11Log = zeros(prams.ntime,ntra);
F12Log = zeros(prams.ntime,ntra);
F21Log = zeros(prams.ntime,ntra);
F22Log = zeros(prams.ntime,ntra);
% allocate memory for positions and deformation gradient

ttLinear = linspace(0,prams.T,prams.ntime);
ttLog = exp(linspace(log(1e-4),log(prams.T),prams.ntime));
%tt = [ttLinear ttLog];
%tt = unique(sort(tt));
%
%sLinear = zeros(prams.ntime,1);
%sLog = zeros(prams.ntime,1);
%
%for k = 1:prams.ntime
%  [~,sLinear(k)] = min(abs(tt - ttLinear(k)));
%  [~,sLog(k)] = min(abs(tt - ttLog(k)));
%end


for k = 1:numel(X0)/2
  message = ['\ntracers ' num2str(2*k/numel(X0)*100,'%04.1f\n') ' %% completed\n'];
  fprintf(message);
  x0 = X0(k);
  y0 = X0(k+numel(X0)/2);
  % initial condition of the kth tracer

  [ttLinear,z] = ode45(odeFun,ttLinear,[x0 y0 1 0 0 1],opts);
  % Solve for position and deformation gradient with ode45 on a linear
  % time scale

  xtraLinear(:,k) = z(:,1);
  ytraLinear(:,k) = z(:,2);
  F11Linear(:,k) = z(:,3);
  F12Linear(:,k) = z(:,4);
  F21Linear(:,k) = z(:,5);
  F22Linear(:,k) = z(:,6);

%  [ttLog,z] = ode45(odeFun,ttLog,[x0 y0 1 0 0 1],opts);
  % Solve for position and deformation gradient with ode45 on a
  % logarithmic time scale
  xtraLog(:,k) = z(:,1);
  ytraLog(:,k) = z(:,2);
  F11Log(:,k) = z(:,3);
  F12Log(:,k) = z(:,4);
  F21Log(:,k) = z(:,5);
  F22Log(:,k) = z(:,6);

  iLeft = true;
  xshifts = [];
  yshifts = [];
  indshifts = [];
  while (iLeft || iRight || iTop || iBottom)
    iLeft = false; iRight = false; iTop = false; iBottom = false;
    if xtraLinear(end,k) > options.xpThresh
      iRight = true;
    elseif xtraLinear(end,k) < options.xmThresh
      iLeft = true;
    elseif ytraLinear(end,k) > options.ypThresh
      iTop = true;
    elseif ytraLinear(end,k) < options.ymThresh
      iBottom = true;
    end

    if (iLeft || iRight || iTop || iBottom)
      if iRight
        indOut = find(xtraLinear(:,k) > options.xpThresh);
      end
      if iLeft
        indOut = find(xtraLinear(:,k) < options.xmThresh);
      end
      if iBottom
        indOut = find(ytraLinear(:,k) < options.ymThresh);
      end
      if iTop
        indOut = find(ytraLinear(:,k) > options.ypThresh);
      end

      indOut = max(max(indOut(1)) - 1,1);
      vel = odeFun(0,[xtraLinear(indOut,k) ytraLinear(indOut,k) ...
          0 0 0 0]);

      l2Dist = (vel(1) - uWin).^2 + (vel(2) - vWin).^2;
      [~,smin] = min(l2Dist(l2Dist>1e-5));
      % don't want to get the same point twice.  This can happen if a
      % tracer goes out the top or bottom of the geometry in the first
      % small section where the tracer is shited in the periodic 'hack' 
      
      z0 = [xWin(smin) yWin(smin) ...
          F11Linear(indOut,k) F12Linear(indOut,k) ...
          F21Linear(indOut,k) F22Linear(indOut,k)];

      indshifts = [indshifts indOut];
      xshifts = [xshifts xtraLinear(indOut,k) - xWin(smin)];
      yshifts = [yshifts ytraLinear(indOut,k) - yWin(smin)];

      [~,z] = ode45(odeFun,ttLinear(indOut:end),...
          z0,opts); 
      if numel(ttLinear) == indOut + 1
        z = [z(1,:) z(end,:)];
      end
      % this happens when indOut is the second to last point of the time
      % points of interest.  If only a 2x1 vector is passed into ode45,
      % it returns the solution at a set of points intermediate to these
      % two points.  We only care about the solution at the first and
      % last point

      xtraLinear(indOut:end,k) = z(:,1);
      ytraLinear(indOut:end,k) = z(:,2);
      F11Linear(indOut:end,k) = z(:,3);
      F12Linear(indOut:end,k) = z(:,4);
      F21Linear(indOut:end,k) = z(:,5);
      F22Linear(indOut:end,k) = z(:,6);
    else
      iLeft = false; iRight = false; iTop = false; iBottom = false;
    end

  end


  if any((xtraLinear(end,k) - centers(:,1)).^2 + ...
      (ytraLinear(end,k) - centers(:,2)).^2 < radii.^2)
    disp('save this index so that I can remove it from the statistics')
    pause
  end


  for j = 1:numel(indshifts)
    xtraLinear(indshifts(j):end,k) = ...
        xtraLinear(indshifts(j):end,k) + xshifts(j);
    ytraLinear(indshifts(j):end,k) = ...
        ytraLinear(indshifts(j):end,k) + yshifts(j);
  end

  if mod(k,1) == 1
    om.writeTracerPositions(ttLinear,xtraLinear(:,1:k),ytraLinear(:,1:k),'linear');
%    om.writeTracerPositions(ttLog,xtraLog(:,1:k),ytraLog(:,1:k),'log');
    if options.defGradient
      om.writeDeformationGradient(time,F11(:,1:k),F12(:,1:k),...
          F21(:,1:k),F22(:,1:k));
    end
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

om.writeTracerPositions(ttLinear,xtraLinear,ytraLinear,'linear');
%om.writeTracerPositions(ttLog,xtraLog,ytraLog,'log');
% one final save at the very end



