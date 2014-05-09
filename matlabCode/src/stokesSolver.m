function stokesSolver(Xinner,Xouter,options,prams)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% DESCRIPTION:
% tracers = stokesSolver(Xouter,Xinner,options,prams) returns the time
% history of a set of tracers
% INPUTS
% Xouter - parameterization of the outer boundary
% Xinner - parameterization of the inner boundaries
% prams - set of parameters
% options - set of options
% OUTPUTS:
% none 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
om = monitor(options,prams);
% object for doing I/O
om.welcomeMessage(options,prams);
% output parameters
om.clearFiles;
% empty the dat and bin files

innerGeom = capsules(Xinner,'inner');
outerGeom = capsules(Xouter,'outer');
clear Xinner
clear Xouter
% create objects for the inner and outer boundaries.
% The outer boundary will need more points so need
% two classes.  Also, think that the double-layer potential
% may be more appropriate for the outer boundary as the
% circle may not be a very good preconditioner.  This way, we can use
% multigrid with Picard

%rhs = [innerGeom.u(:); outerGeom.u(:)];
rhs = [zeros(2*innerGeom.N*innerGeom.nv,1); outerGeom.u(:)];
% right-hand side which corresponds to no-slip on the solid walls

op = poten(prams.Ninner,options.fmm);
% object for evaluating layer potentials

tic
[eta,flag,relres,iter,relresvec] = gmres(...
    @(X) op.matVecMultiply(X,innerGeom,outerGeom),...
    rhs,[],prams.gmresTol,prams.maxIter,...
    @(X) op.matVecPreco(X,innerGeom));
%[eta,flag,relres,iter,relresvec] = gmres(...
%    @(X) op.matVecMultiply(X,innerGeom,outerGeom),...
%    rhs,[],prams.gmresTol,prams.maxIter);
% do unpreconditioned GMRES to find density function
om.writeStars
message = ['****     GMRES took ' num2str(toc,'%4.2e') ...
    ' seconds     ****'];
om.writeMessage(message,'%s\n');
if (flag ~= 0)
  message = 'GMRES HAD A PROBLEM';
else
  message = ['****     GMRES required ' num2str(iter(2),'%3d'),...
      ' iterations    ****'];
end
om.writeMessage(message,'%s\n');
om.writeStars
om.writeMessage(' ');

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

om.writeGeometry(innerGeom.X,outerGeom.X,sigmaInner,sigmaOuter);
% write the density function to the data file for post processing
om.writeStars
message = '****    Density function and geometry    ****';
om.writeMessage(message);
message = '****         written to bin file         ****';
om.writeMessage(message);
om.writeStars
om.writeMessage(' ');




