addpath ../src

load radii2.dat;
%load centers2.dat;
load centers3.dat;

radii2 = radii2(1:10);
%centers2 = centers2(1:10,:);
centers3 = centers3(1:10,:);

prams.Nouter = 1024;
% number of points on outer solid wall
prams.Ninner = 256;
% number of points per circle exclusion
prams.nv = numel(radii2);
% number of exclusions
prams.gmresTol = 1e-8;
% gmres tolerance
prams.maxIter = min(2*(prams.Nouter + prams.nv*prams.Ninner),1500);
% maximum number of gmres iterations
%prams.atol = 1e-6;
%prams.rtol = 1e-3;
prams.atol = 1e-9;
prams.rtol = 1e-6;
% absolute and relative tolerances for ode45
%prams.T = 30*3;
%% time horizon for ode45
%prams.ntime = 1500*3 + 1;
%% number of time steps that ode45 will output
prams.T = 1e0;
prams.ntime = 100;

% Different options
options.bieSolve = false; 
options.computeEuler = true;
options.tracersSimulation = true;
options.axis = [-8 38 -0.1 5.3];
options.dataFile = 'output/circles2Data.bin';
options.farField = 'circles';
options.fmm = true;
options.logFile = 'output/circles2.log';
options.profile = false;
options.saveData = true;
options.usePlot = false;
options.verbose = true;

oc = curve;
Xouter = oc.initConfig(prams.Nouter,'square2');
% outer most boundary
%Xinner = oc.initConfig(prams.Ninner,'circles2', ...
%          'nv',prams.nv, ...
%          'center',centers2, ...
%          'radii',radii2);
Xinner = oc.initConfig(prams.Ninner,'circles', ...
          'nv',prams.nv, ...
          'center',centers3, ...
          'radii',radii2);
% built centers3 from centers2 where I did the shifting flipping of the
% centers rather than the geometry.  Then, can do quick checks for
% determing interior and exterior points when computing Eulerian grid
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
  ntra = 1;
  [xtar,ytar] = initialTracers(radii2,centers3,ntra);
  X0 = [xtar(:);ytar(:)];
  % initial tracer locations
  fileName = 'output/circles2Data.bin';
  % file that has all the necessary density function and geometry stored
  options.xmin = 0;
  options.xmax = 35;
  options.nx = 20;
  % min, max, and number of Euler locations in x direction
  options.ymin = 0.2;
  options.ymax = 5;
  options.ny = 20;
  % min, max, and number of Euler locations in y direction
  options.nparts = 1;
  % need to compute in sections otherwise seem to run out of memory
  options.ymThresh = options.ymin + 3;
  options.ypThresh = options.ymax - 2;
  % thresholds where velocity will be set to zero

  [time,xtra,ytra,F11,F12,F21,F22] = tracers(...
      X0,options,prams,fileName);
  % simulate tracers. Each column represents a tracer and each row
  % represents the time variable
end


if options.profile
  profile viewer;
  profile off;
  filename = [options.logFile(1:end-4) 'Profile'];
  profsave(profile('info'),filename);
end
% save the profile


