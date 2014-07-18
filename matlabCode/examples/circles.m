addpath ../src

load radii.dat;
load centers.dat;
%centers = centers(1:50,:);
%radii = radii(1:50);

prams.Nouter = 1024;
% number of points on outer solid wall
prams.Ninner = 256;
% number of points per circle exclusion
prams.nv = numel(radii);
%prams.nv = 1;
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
options.computeEuler = false;
options.tracersSimulation = true;
options.axis = [-0.5 5.5 0 30];
%options.axis = [3.65 3.8 23.95 24.15];
%options.axis = [-0.5 5.5 15 30];
%options.axis = [3.5 4.5 23.8 24.5];
options.axis = [-0.1 5.1 0 30];
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
%  [xtar,ytar] = meshgrid(linspace(0.2,4.8,200),linspace(30,30,1));
%  [xtar,ytar] = meshgrid(linspace(1,4,2),linspace(30,30,1));
%  [xtar,ytar] = meshgrid(linspace(2.3,2.3,1),linspace(30,30,1));
%  xtar = 2.0; ytar = 30;
  ntra = 1000;
  [xtar,ytar] = initialTracers(radii,centers,ntra);
  X0 = [xtar(:);ytar(:)];
%  X0 = [1.7;2.6;22;18.4];
  % initial tracer locations
  fileName = 'output/circlesData.bin';
  % file that has all the necessary density function and geometry stored
  options.xmin = 0.25;
  options.xmax = 4.53;
  options.nx = 800;
%  options.xmin = 1.6;
%  options.xmax = 1.8;
%  options.nx = 201;
%  dx = (1.25e-2)/8;
%  options.xmin = 1.7 - dx;
%  options.xmax = 1.7 + dx;
%  options.nx = 11;
  % min, max, and number of Euler locations in x direction
  options.ymin = 0;
  options.ymax = 33;
  options.ny = 7200;
  options.nparts = 5;
%  options.ymin = 21.7;
%  options.ymax = 21.9;
%  options.ny = 201;
%  options.ymin = 22 - dx;
%  options.ymax = 22 + dx;
%  options.ny = 11;
  % min, max, and number of Euler locations in y direction
%  options.nparts = 1;
  % need to compute in sections otherwise seem to run out of memory
  options.ymThresh = options.ymin + 3;
  options.ypThresh = options.ymax - 2;
%  options.ymThresh = -100;
%  options.ypThresh = 100;
  % thresholds where velocity will be set to zero

  [time,xtra,ytra,F11,F12,F21,F22] = tracers(...
      X0,options,prams,fileName);
  % simulate tracers. Each column represents a tracer and each row
  % represents the time variable
end


if options.profile
  profile off;
  filename = [options.logFile(1:end-4) 'Profile'];
  profsave(profile('info'),filename);
end
% save the profile


