addpath ../src

load radii.dat;
load centers.dat;

prams.Nouter = 512;
% number of points on outer solid wall
prams.Ninner = 128;
%prams.Ninner = prams.Nouter;
% number of points per circle exlcusions
%prams.nv = numel(radii);
prams.nv = 150;
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

options.farField = 'circles';
options.fmm = true;
options.verbose = true;
options.dataFile = 'output/circlesData.bin';
options.logFile = 'output/circles.log';
options.usePlot = false;
options.profile = false;

oc = curve;
Xouter = oc.initConfig(prams.Nouter,'square');
% outer most boundary
Xinner = oc.initConfig(prams.Ninner,'circles', ...
          'nv',prams.nv, ...
          'center',centers, ...
          'radii',radii);
% circular exclusions

%[xtar,ytar] = meshgrid(linspace(0.1,4.9,40),linspace(-6,35,80));
%[xtar,ytar] = meshgrid(linspace(3.5,4.5,50),linspace(23.5,24.5,50));
[xtar,ytar] = meshgrid(linspace(0.5,4.5,10),linspace(35,35,1));
X0 = [xtar(:);ytar(:)];
% initial conditions
prams.Ntracers = numel(xtar);

if options.profile
  profile off; profile on;
end

[time,tracers] = stokesSolver(X0,Xouter,Xinner,prams,options);

if options.profile
  profile off;
  filename = [options.logFile(1:end-4) 'Profile'];
  profsave(profile('info'),filename);
end


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


