%clear all
addpath ../src

load radii12.dat;
load centers12.dat;
radii = radii12;
centers = centers12;

%radii = radii(1:300);
%centers = centers(1:300,:);

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
prams.atol = 1e-9;
prams.rtol = 1e-6;
%prams.atol = 1e-6;
%prams.rtol = 1e-3;
% absolute and relative tolerances for ode45
%prams.T = 30*3;
%% time horizon for ode45
%prams.ntime = 1500*3 + 1;
%% number of time steps that ode45 will output
prams.T = 3.2306e1;
prams.ntime = 200;
%prams.T = 40;
%prams.ntime = 4000;

% Different options
options.bieSolve = false; 
options.computeEuler = false;
options.tracersSimulation = true;
options.defGradient = false;
options.axis = [-8.0 44 -0.2 5.4];
options.dataFile = '/scratch/quaife/porousSimulations/results/newGeoms/circles12Data.bin';
%options.dataFile = 'output/circles12Data.bin';
options.farField = 'circles';
options.fmm = false;
options.logFile = 'output/circles12.log';
options.profile = false;
options.saveData = true;
options.usePlot = false;
options.verbose = true;

oc = curve;
Xouter = oc.initConfig(prams.Nouter,'square12');
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
  options.xmThresh = 1;
  options.xpThresh = 37;
  options.ymThresh = 0.5;
  options.ypThresh = 4.7;
%  options.ymThresh = 2.0;
%  options.ypThresh = 3.0;
  % thresholds where velocity will be set to zero

  ntra = 1;
%  [xtar,ytar] = initialTracers(radii,centers,ntra,options);
%  X0 = [xtar(:);ytar(:)];
  X0 = [6.610844675680466e0;4.665988941855113e0];
%  X0 = [6.331921348743297e0;4.398071171843299e0];
%  X0 = [3.461225119059567e1;2.713024928665051e+00];
  % initial tracer locations
  fileName = options.dataFile;
  % file that has all the necessary density function and geometry stored
  options.xmin = -1;
  options.xmax = 39;
  options.nx = 10000;
  % min, max, and number of Euler locations in x direction
  options.ymin = 0.001;
  options.ymax = 5.199;
  options.ny = 1000;
  % min, max, and number of Euler locations in y direction
  options.nparts = 100;
  % need to compute in sections otherwise seem to run out of memory

  tracers(X0,options,prams,radii,centers,fileName);
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


