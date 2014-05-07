addpath ../src

load radii.dat;
load centers.dat;

prams.Nouter = 1024;
% number of points on outer solid wall
prams.Ninner = 256;
% number of points per circle exclusion
prams.nv = 2;
% number of exclusions
prams.gmresTol = 1e-8;
% gmres tolerance
prams.maxIter = 2*(prams.Nouter + prams.nv*prams.Ninner);
%prams.maxIter = 1;
% maximum number of gmres iterations
prams.atol = 1e-6;
prams.rtol = 1e-3;
% absolute and relative tolerances for ode45
prams.T = 150;
% time horizon for ode45
prams.ntime = 1501;
% number of time steps that ode45 will output

% Different options
options.bieSolve = true;
options.computeEuler = true;
options.tracersSimulation = false;
options.axis = [-0.5 5.5 0 30];
%options.axis = [3.5 4.5 23.8 24.5];
%options.axis = [-0.1 5.1 20 30];
options.dataFile = 'output/circlesData.bin';
options.farField = 'circles';
options.fmm = true;
options.logFile = 'output/circles.log';
options.profile = false;
options.saveData = true;
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
  [xtar,ytar] = meshgrid(linspace(0.2,4.8,200),linspace(30,30,1));
%  [xtar,ytar] = meshgrid(linspace(0.2,4.8,10),linspace(30,30,1));
  X0 = [xtar(:);ytar(:)];
  % initial tracer locations
  fileName = 'output/circlesData.bin';
  % file that has all the necessary density function and geometry stored
  options.xmin = 0.05;
  options.xmax = 4.95;
  options.nx = 200;
  % min, max, and number of Euler locations in x direction
  options.ymin = 0;
  options.ymax = 35;
  options.ny = 1800;
  % min, max, and number of Euler locations in y direction
  options.ymThresh = options.ymin + 2;
  options.ypThresh = options.ymax - 2;
%  options.ymThresh = 23.6;
%  options.ypThresh = 24.4;
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


