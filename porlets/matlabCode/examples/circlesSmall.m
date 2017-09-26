%clear all
addpath ../src

%radii = [0.5 0.2 0.5 0.4];
%centers = [[2 2];[4 3];[3 1];[5 2]];
%radii = [0.1];
%centers = [15 2.5];
load radii28.dat;
load centers28.dat;

radii = radii28(1:20);
centers = centers28(1:20,:);

prams.Nouter = 2048;
% number of points on outer solid wall
prams.Ninner = 256;
% number of points per circle exclusion
prams.nv = numel(radii);
% number of exclusions
prams.gmresTol = 1e-6;
% gmres tolerance
prams.maxIter = min(2*(prams.Nouter + prams.nv*prams.Ninner),5000);
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
prams.T = 10;
prams.ntime = 150;

% Different options
options.bieSolve = true; 
options.computeEuler = true;
options.tracersSimulation = true;
options.defGradient = false;
options.axis = [-6.3 36.7 -0.2 5.4];
options.dataFile = 'output/circlesSmallData.bin';
options.farField = 'circles';
options.fmm = false;
options.logFile = 'output/circlesSmall.log';
options.profile = false;
options.saveData = true;
options.usePlot = false;
options.verbose = true;

oc = curve;
Xouter = oc.initConfig(prams.Nouter,'square28');
% outer most boundary
Xinner = oc.initConfig(prams.Ninner,'circles', ...
          'nv',prams.nv, ...
          'center',centers, ...
          'radii',radii);
% built centers3 from centers2 where I did the shifting flipping of the
% centers rather than the geometry.  Then, can do quick checks for
% determing interior and exterior points when computing Eulerian grid
% circular exclusions

%figure(1); clf; hold on
%plot(Xouter(1:end/2),Xouter(end/2+1:end),'k')
%axis equal;
%fill(Xinner(1:end/2,:),Xinner(end/2+1:end,:),'k');
%axis(options.axis)
%pause

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
  ntra = 0;
  [xtar,ytar] = initialTracers(radii,centers,ntra,options);
  X0 = [xtar(:);ytar(:)];
%  X0 = [];
% X0 = [30;4.5];
  % initial tracer locations
  fileName = options.dataFile;
  % file that has all the necessary density function and geometry stored
  options.xmin = 0;
  options.xmax = 30;
  options.nx = 200;
  % min, max, and number of Euler locations in x direction
  options.ymin = 0.1;
  options.ymax = 5.1;
  options.ny = 200;
  % min, max, and number of Euler locations in y direction
  options.nparts = 5;
  % need to compute in sections otherwise seem to run out of memory
  options.xmThresh = options.xmin + 0;
  options.xpThresh = options.xmax - 0;
  options.ymThresh = options.ymin + 0;
  options.ypThresh = options.ymax - 0;
  % thresholds where velocity will be set to zero

  tracers(X0,options,prams,radii,centers,fileName)
%  [time,xtra,ytra,F11,F12,F21,F22] = tracers(...
%      X0,options,prams,fileName);
  % simulate tracers. Each column represents a tracer and each row
  % represents the time variable
end


if options.profile
%  profile viewer;
  profile off;
  filename = [options.logFile(1:end-4) 'Profile'];
  profsave(profile('info'),filename);
end
% save the profile


