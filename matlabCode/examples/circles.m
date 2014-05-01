addpath ../src

load radii.dat;
load centers.dat;

prams.Nouter = 512;
% number of points on outer solid wall
prams.Ninner = 128;
% number of points per circle exclusion
prams.nv = 1;
% number of exclusions
prams.gmresTol = 1e-8;
% gmres tolerance
prams.maxIter = 2*(prams.Nouter + prams.nv*prams.Ninner);
%prams.maxIter = 1;
% maximum number of gmres iterations
prams.atol = 1e-6;
prams.rtol = 1e-3;
% absolute and relative tolerances for ode45
prams.T = 50;
% time horizon for ode45
prams.ntime = 101;
% number of time steps that ode45 will output

% Different options
options.bieSolve = false;
options.dataFile = 'output/circlesData.bin';
options.farField = 'circles';
options.fmm = false;
options.logFile = 'output/circles.log';
options.profile = false;
options.saveData = true;
options.tracersSimulation = true;
options.usePlot = false;
options.verbose = true;

oc = curve;
Xouter = oc.initConfig(prams.Nouter,'square');
% outer most boundary
Xinner = oc.initConfig(prams.Ninner,'circles', ...
          'nv',prams.nv, ...
          'center',centers, ...
          'radii',radii);
% circular exclusions

if options.profile
  profile off; profile on;
end

if options.bieSolve
  stokesSolver(Xinner,Xouter,options,prams);
end
% solve density function and write to .bin files.  It this calculation
% has already ben done, everything can be loaded in tracers to do the
% Lagrange tracker simulation.




if options.tracersSimulation
%  [xtar,ytar] = meshgrid(linspace(0.5,4.5,100),linspace(30,35,2));
  [xtar,ytar] = meshgrid(linspace(3.7,4.1,100),linspace(24.3,24.3,1));
  X0 = [xtar(:);ytar(:)];
  % initial tracer locations
  fileName = 'output/circlesData.bin';
  % file that has all the necessary density function and geometry stored
  options.xmin = 3.6; %0.1;
  options.xmax = 4.2;
  options.nx = 200;
  % min, max, and number of Euler locations in x direction
  options.ymin = 23.8; %-10;
  options.ymax = 24.4; %35;
  options.ny = 200; %1000;
  % min, max, and number of Euler locations in y direction
  options.ymThresh = 0;
  options.ypThresh = 40;
  % thresholds where velocity will be set to zero

  [time,xtra,ytra] = tracers(X0,options,prams,fileName);
  % simulate tracers. Each column represents a tracer and each row
  % represents the time variable
end



if options.profile
  profile off;
  filename = [options.logFile(1:end-4) 'Profile'];
  profsave(profile('info'),filename);
end
% save the profile


%if options.usePlot
%  figure(1); clf;
%  hold on
%  plot(Xouter(1:end/2,:),Xouter(end/2+1:end,:),'k')
%  axis equal
%
%  oc = curve;
%  [tracerx,tracery] = oc.getXY(tracers');
%  plot(tracerx',tracery','b-')
%
%  fill(Xinner(1:end/2,:),Xinner(end/2+1:end,:),'k')
%end


