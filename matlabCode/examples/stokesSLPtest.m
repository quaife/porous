%clear all
addpath ../src

load radii.dat;
load centers.dat;

prams.Nouter = [];
prams.Ninner = 256;
% number of points per circle exclusion
prams.nv = 465;
% number of exclusions
prams.gmresTol = 1e-8;
% gmres tolerance
prams.maxIter = min(2*prams.nv*prams.Ninner,500);
%prams.maxIter = 1;
%prams.maxIter = 30;
% maximum number of gmres iterations

% Different options
options.bieSolve = true;
options.farField = 'circles';
options.fmm = true;
options.profile = false;
options.saveData = true;
options.verbose = true;
options.dataFile = ' ';
options.logFile = 'output/stokesSLPtest.log';
options.axis = [];

theta = (0:prams.Ninner-1)'*2*pi/prams.Ninner;
%X = [cos(theta) cos(theta)+2.01 cos(theta) cos(theta)+2.01; ...
%     sin(theta) sin(theta) sin(theta)+2.01 sin(theta)+2.01];
%Xinner = 2*[cos(theta);sin(theta)];
%Xinner = [cos(theta) cos(theta)+2.1; ...
%     sin(theta) sin(theta)];
oc = curve;
Xinner = oc.initConfig(prams.Ninner,'circles', ...
          'nv',prams.nv, ...
          'center',centers, ...
          'radii',radii);
% circular exclusions
XinnerCoarse = oc.initConfig(32,'circles',...
          'nv',prams.nv, ...
          'center',centers, ...
          'radii',radii);
% coarse grid

%if options.profile
%  profile off; profile on;
%end

[S,P] = SLPsolver(Xinner,XinnerCoarse,options,prams);
% solve density function and write to .bin files.  It this calculation
% has already ben done, everything can be loaded in tracers to do the
% Lagrange tracker simulation.



%if options.profile
%  profile off;
%  filename = [options.logFile(1:end-4) 'Profile'];
%  profsave(profile('info'),filename);
%end
%% save the profile


