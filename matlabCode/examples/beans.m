clear all
addpath ../src

load radii.dat;
load centers.dat;
load radiiBeans.dat;
load centersBeans.dat;

prams.Nouter = 2048/16;
% number of points on outer solid wall
prams.Ninner = 256/4;
% number of points per circle exclusion
prams.nv = numel(radii);
% number of exclusions
prams.gmresTol = 1e-6;
% gmres tolerance
prams.maxIter = min(2*(prams.Nouter + prams.nv*prams.Ninner),500);
% maximum number of gmres iterations
prams.atol = 1e-9;
prams.rtol = 1e-6;
% absolute and relative tolerances for ode45
prams.T = 1e0;
% time horizon for ode45
prams.ntime = 220;
% number of time steps that ode45 will output

% Different options
options.bieSolve = true; 
options.computeEuler = false;
options.tracersSimulation = false;
options.defGradient = false;
options.axis = [-8 38 -0.1 5.3];
options.dataFile = 'output/BeansData.bin';
options.farField = 'circles';
options.fmm = true;
options.logFile = 'output/beans.log';
options.profile = false;
options.saveData = true;
options.usePlot = true;
options.verbose = true;

oc = curve;
Xouter = oc.initConfig(prams.Nouter,'square');
% outer most boundary
Xinner = oc.initConfig(prams.Ninner,'beans', ...
          'nv',prams.nv, ...
          'center',centers, ...
          'radii',radii,...
          'centerBeans',centersBeans, ...
          'radiiBeans',radiiBeans);
% built centers and centersBenas where I did the shifting flipping of
% the centers rather than the geometry.  Then, can do quick checks for
% determing interior and exterior points when computing Eulerian grid
% circular exclusions
Xinner = [Xinner(:,326:326)];
prams.nv = size(Xinner,2);
size(Xinner)

%figure(2); clf; hold on
%plot(Xouter(1:end/2),Xouter(end/2+1:end),'k')
%axis equal;
%fill(Xinner(1:end/2,340:465),Xinner(end/2+1:end,340:465),'k');
%%fill(Xinner(1:end/2,466:end),Xinner(end/2+1:end,466:end),'w','edgecolor','w');
%fill(Xinner(1:end/2,1:339),Xinner(end/2+1:end,1:339),'k');




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
  ntra = 10000;
  [xtar,ytar] = initialTracers(radii,centers,ntra);
  X0 = [xtar(:);ytar(:)];
%  X0 = [];
% X0 = [30;4.5];
  % initial tracer locations
  fileName = options.dataFile;
  % file that has all the necessary density function and geometry stored
  options.xmin = 0;
  options.xmax = 35;
  options.nx = 9000;
  % min, max, and number of Euler locations in x direction
  options.ymin = 0.001;
  options.ymax = 5.199;
  options.ny = 1000;
  % min, max, and number of Euler locations in y direction
  options.nparts = 10;
  % need to compute in sections otherwise seem to run out of memory
  options.xmThresh = options.xmin + 0;
  options.xpThresh = options.xmax - 0;
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


