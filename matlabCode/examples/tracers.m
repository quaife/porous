addpath ../src

[Ninner,Nouter,nv,Xinner,Xouter,sigmaInner,sigmaOuter] = ...
    loadFile('output/circlesData.bin');
% load all the information about the geometry and density function
innerGeom = capsules(Xinner,'inner');
outerGeom = capsules(Xouter,'inner');
% build objects for the inner and outer boundaries

fmm = true;
op = poten(Ninner,fmm);

%[eulerX,eulerY] = meshgrid(0.1:0.4:4.9,0:0.4:32);
[eulerX,eulerY] = meshgrid(linspace(0.1,4.9,100),linspace(-10,40,1000));

% choose an Eulerian grid

vel = op.layerEval(0,[eulerX(:);eulerY(:)],...
    innerGeom,outerGeom,sigmaInner,sigmaOuter);
% evalute velocity on an Eulerian grid

u = reshape(vel(1:end/2),size(eulerX));
v = reshape(vel(end/2+1:end),size(eulerY));
% put velocity field in format that works well for interp2

[xtra,ytra] = meshgrid(linspace(0.4,4.6,1000),[30]);
xtra = xtra(:); ytra = ytra(:);
X0 = [xtra;ytra];
% initial condition of tracers


odeFun = @(t,z) op.interpolateLayerPot(t,z,eulerX,eulerY,u,v);
[t,Xtra] = ode45(odeFun,linspace(0,20,2000),X0);

