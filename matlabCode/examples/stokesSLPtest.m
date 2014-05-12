%clear all
addpath ../src

load radii.dat;
load centers.dat;

prams.Nouter = [];
prams.Ninner = 8;
% number of points per circle exclusion
prams.nv = 100;
% number of exclusions
prams.gmresTol = 1e-8;
% gmres tolerance
prams.maxIter = min(2*prams.nv*prams.Ninner,200);
% maximum number of gmres iterations

% Different options
options.bieSolve = true;
options.farField = 'circles';
options.fmm = false;
options.profile = false;
options.saveData = true;
options.verbose = true;
options.dataFile = ' ';
options.logFile = 'output/stokesSLPtest.log';
options.axis = [];

theta = (0:prams.Ninner-1)'*2*pi/prams.Ninner;
%X = [1/2*cos(theta) cos(theta)+3 cos(theta) cos(theta)+3; ...
%     1/2*sin(theta) sin(theta) sin(theta)+3 sin(theta)+3];
%Xinner = 2*[cos(theta);sin(theta)];
%Xinner = [cos(theta) cos(theta)+3; ...
%     sin(theta) sin(theta)];
oc = curve;
Xinner = oc.initConfig(prams.Ninner,'circles', ...
          'nv',prams.nv, ...
          'center',centers, ...
          'radii',radii);
% circular exclusions

if options.profile
  profile off; profile on;
end

[S,P] = SLPsolver(Xinner,options,prams);
% solve density function and write to .bin files.  It this calculation
% has already ben done, everything can be loaded in tracers to do the
% Lagrange tracker simulation.



if options.profile
  profile off;
  filename = [options.logFile(1:end-4) 'Profile'];
  profsave(profile('info'),filename);
end
% save the profile


