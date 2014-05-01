function [time,tracers] = stokesSolver(X0,Xouter,Xinner,prams,options)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% DESCRIPTION:
% tracers = stokesSolver(Xouter,Xinner,prams,options) returns the time
% history of a set of tracers
% INPUTS:
% tracers - time history of a set of tracers
% OUTPUTS
% X0 - initial conditions for tracers
% Xouter - parameterization of the outer boundary
% Xinner - parameterization of the inner boundaries
% prams - set of parameters
% options - set of options
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
om = monitor(options,prams);
% object for doing I/O
om.welcomeMessage(options,prams);
% output parameters

innerGeom = capsules(Xinner,'inner');
outerGeom = capsules(Xouter,'outer');
% create objects for the inner and outer boundaries.
% The outer boundary will need more points so need
% two classes.  Also, think that the double-layer potential
% may be more appropriate for the outer boundary as the
% circle may not be a very good preconditioner.  This way, we can use
% multigrid with Picard

rhs = [innerGeom.u(:); 1*outerGeom.u(:)];
% right-hand side which corresponds to no-slip on the solid walls

op = poten(prams.Ninner,options.fmm);
% object for evaluating layer potentials

tic
[eta,flag,relres,iter,relresvec] = gmres(...
    @(X) op.matVecMultiply(X,innerGeom,outerGeom),...
    rhs,[],prams.gmresTol,prams.maxIter,...
    @(X) op.matVecPreco(X,innerGeom));
% do preconditioned GMRES to find density function
om.writeStars
message = ['GMRES took ' num2str(toc,'%4.2e') ' seconds'];
om.writeMessage(message,'%s\n');
if (flag ~= 0)
  message = 'GMRES HAD A PROBLEM';
else
  message = ['GMRES required ' num2str(iter(2),'%3d'),' iterations'];
end
om.writeMessage(message,'%s\n');
om.writeStars


sigmaInner = zeros(2*prams.Ninner,prams.nv);
% initialize space for desnity function along inner boundaries
for k = 1:prams.nv
  istart = 2*(k-1)*prams.Ninner + 1;
  iend = istart + 2*prams.Ninner - 1;
  sigmaInner(:,k) = eta(istart:iend);
end
% unstack the density function at the inner boundaries

istart = 2*prams.nv*prams.Ninner + 1;
iend = istart + 2*prams.Nouter - 1;
sigmaOuter = eta(istart:iend);
% unstack the density function at the outer boundary

om.writeOutput(Xinner,Xouter,sigmaInner,sigmaOuter);
% write the density function to the data file for post processing


timeIntegrator = false;
if timeIntegrator
  odeFun = @(t,z) ...
      op.layerEval(t,z,innerGeom,outerGeom,sigmaInner,sigmaOuter);
  %vel = odeFun(0,X0);
  % compute the velocity at the intial points for forming a quiver plot
  opts.AbsTol = prams.rtol;
  opts.RelTol = prams.atol;
  % tolerances for ODE solver
  T = prams.T; 
  ntime = prams.ntime;
  % time horizon

  tic
  [time,tracers] = ode45(odeFun,linspace(0,T,ntime),X0,opts);
  om.writeStars
  message = ['ode45 took ' num2str(toc,'%4.2e') ' iterations'];
  om.writeMessage(message,'%s\n');
else
  time = [];
  tracers = [];
end
% find tracer positions
