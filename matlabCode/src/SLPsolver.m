function [S,P] = SLPsolve(Xinner,XinnerCoarse,options,prams)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% DESCRIPTION:
% INPUTS
% Xinner - parameterization of the inner boundaries
% XinnerCoarse - parameterization of the inner boundary on a coarser
% grid
% options - set of options
% prams - set of parameters
% OUTPUTS:
% S and P are the matrix and the preconditioner which can be formed for
% suffuciently small problems
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
global matVecLarge matVecSmall
matVecLarge = 0;
matVecSmall = 0;

if 0
  preco = 'none';
  % unpreconditioned gmres
end
if 0
  preco = 'BD';
  % block-diagonal preconditioned gmres
end
if 0
  preco = '2grid';
  % V(1,0) preconditioned gmres
end
if 1
  preco = '2gridIter';
  % V(1,0) iterative (no gmres) 
end

op = poten(prams.Ninner,options.fmm);
% object for evaluating layer potentials
om = monitor(options,prams);
innerGeom = capsules(Xinner,'inner');
NearI2I = innerGeom.getZone(Xinner,1);
if strcmp(preco(1:2),'2g')
  innerGeomCoarse = capsules(XinnerCoarse,'inner');
  NearI2ICoarse = innerGeomCoarse.getZone(XinnerCoarse,1);
end
% create objects for the inner and outer boundaries.  Also create the
% near-singular integration structures

theta = (0:prams.Ninner-1)'*2*pi/prams.Ninner;
etaTrue = zeros(2*prams.Ninner,prams.nv);
for k = 1:prams.nv
  etaTrue(:,k) = [exp(cos(theta));exp(sin(theta))];
end
rhs = op.SLPmatVecMultiply(etaTrue,innerGeom,NearI2I);
rhs = reshape(rhs,2*prams.Ninner,prams.nv);
% right-hand side which corresponds to no-slip on the solid walls

matVecLarge = 0;
matVecSmall = 0;
% reset to 0 since we are now doing solver

if options.profile
  profile off; profile on;
end

tic
if strcmp(preco,'BD')
  fprintf('Using block-diagonal preconditioner\n')
  [eta1,flag,relres,iter,relresvec1] = gmres(...
      @(X) op.SLPmatVecMultiply(X,innerGeom,NearI2I),...
      rhs(:),[],prams.gmresTol,prams.maxIter,...
      @(X) op.matVecInvBD(X,innerGeom,[]));

elseif strcmp(preco,'2grid')
  fprintf('Using two grid V-cycle\n')
  maxIter = min(prams.maxIter,2*innerGeomCoarse.N*innerGeom.nv);
  % maximum number of iterations done at the coarse grid solve
  [eta1,flag,relres,iter,relresvec1] = gmres(...
      @(X) op.SLPmatVecMultiply(X,innerGeom,NearI2I),...
      rhs(:),[],prams.gmresTol,prams.maxIter,...
      @(X) op.twoGridVcycle(X,innerGeom,innerGeomCoarse,...
        NearI2I,NearI2ICoarse,1e-10,maxIter));

elseif strcmp(preco,'2gridIter')
  eta1 = rhs;
  normRes = 1e10;
  iter = [0 0];
  % keep track of number of iterations taken using same structure that
  % is used by gmres
  maxVcycles = 100;
  % maximum number of allowable V-cycles
  maxIter = 10;
  while (normRes > prams.gmresTol && iter(2) < maxVcycles)
    iter(2) = iter(2) + 1;
    [eta1,normRes,gmresIter] = op.twoGridVcycle(rhs(:),...
        innerGeom,innerGeomCoarse,NearI2I,NearI2ICoarse,...
        1e-2,maxIter,eta1(:));
    maxIter = max(10,ceil(1.5*gmresIter));
    % set number of allowable coarse grid iterations to be 50% more
    % than the previous number of iterations required, as long as this
    % number exeeds 10
    message = ['Coarse grid solve used ' num2str(iter) ...
        ' of the allowable ' num2str(maxIter) ' iterations'];
    om.writeMessage(message);
    message = ['Iteration ' num2str(iter(2)) ...
        ' residual is ' num2str(normRes,'%4.2e')];
    om.writeMessage(message);
    om.writeStars
    eta1 = reshape(eta1,2*innerGeom.N,innerGeom.nv);
  end
  flag = 0;

else
  fprintf('Not using preconditioner\n')
  [eta1,flag,relres,iter,relresvec1] = gmres(...
      @(X) op.SLPmatVecMultiply(X,innerGeom,NearI2I),...
      rhs(:),[],prams.gmresTol,prams.maxIter);
end
% find density function using gmres, block-diagonal preconditioned
% gmres, V-cycle preconditioned gmres, or V-cycle iterations

eta1 = reshape(eta1,2*innerGeom.N,innerGeom.nv);
fprintf('Number of large matvecs is %d\n',matVecLarge)
fprintf('Number of small matvecs is %d\n',matVecSmall)
res = norm(op.SLPmatVecMultiply(eta1,innerGeom,NearI2I)...
    -rhs(:))/norm(rhs(:));
fprintf('Final relative residual is %4.2e\n',res)
err = norm(eta1 - etaTrue)/norm(etaTrue);
fprintf('Final relative error is %4.2e\n',err)

om.writeStars
message = ['****    pGMRES took ' num2str(toc,'%4.2e') ...
    ' seconds     ****'];
om.writeMessage(message,'%s\n');
if (flag == 0)
  message = ['****    pGMRES required ' num2str(iter(2),'%3d'),...
      ' iterations    ****'];
  om.writeMessage(message,'%s\n');
elseif (flag == 1)
  message = '****    GMRES tolerance not achieved     ****';
  om.writeMessage(message,'%s\n');
  message = ['****    achieved tolerance is ' ...
      num2str(relres,'%4.2e') '   ****'];
  om.writeMessage(message,'%s\n');
  message = ['****    pGMRES took ' num2str(iter(2),'%3d'),...
      ' iterations        ****'];
  om.writeMessage(message,'%s\n');
else 
  message = 'GMRES HAD A PROBLEM';
  om.writeMessage(message,'%s\n');
end
om.writeStars
om.writeMessage(' ');

if options.profile
  profile off;
  filename = [options.logFile(1:end-4) 'Profile'];
  profsave(profile('info'),filename);
end
% save the profile

%sa = innerGeom.sa;
%rhs2 = rhs.*sqrt(2*pi*[sa;sa])/sqrt(innerGeom.N);
%
%tic
%if strcmp(preco,'BD')
%  [eta2,flag,relres,iter,relresvec2] = minres(...
%      @(X) op.SLPmatVecMultiplyGalerkin(X,innerGeom),...
%      rhs2(:),prams.gmresTol,prams.maxIter,...
%      @(X) op.matVecInvBD(X,innerGeom,[]));
%else
%  [eta2,flag,relres,iter,relresvec2] = minres(...
%      @(X) op.SLPmatVecMultiplyGalerkin(X,innerGeom),...
%      rhs2(:),prams.gmresTol,prams.maxIter);
%end
%% do MINRES to find density function
%om.writeStars
%message = ['****    minres took ' num2str(toc,'%4.2e') ...
%    ' seconds     ****'];
%om.writeMessage(message,'%s\n');
%if (flag == 0)
%  message = ['****    minres required ' num2str(iter,'%3d'),...
%      ' iterations    ****'];
%  om.writeMessage(message,'%s\n');
%elseif (flag == 1)
%  message = '****    minres tolerance not achieved    ****';
%  om.writeMessage(message,'%s\n');
%  message = ['****    achieved tolerance is ' ...
%      num2str(relres,'%4.2e') '   ****'];
%  om.writeMessage(message,'%s\n');
%  message = ['****    minres took ' num2str(iter,'%3d'),...
%      ' iterations        ****'];
%  om.writeMessage(message,'%s\n');
%else 
%  message = 'MINRES HAD A PROBLEM';
%  om.writeMessage(message,'%s\n');
%  flag
%end
%om.writeStars
%om.writeMessage(' ');
%
%
%tic
%if strcmp(preco,'BD')
%  [eta3,flag,relres,iter,relresvec3] = symmlq(...
%      @(X) op.SLPmatVecMultiplyGalerkin(X,innerGeom),...
%      rhs2(:),prams.gmresTol,prams.maxIter,...
%      @(X) op.matVecInvBD(X,innerGeom,[]));
%else
%  [eta3,flag,relres,iter,relresvec3] = symmlq(...
%      @(X) op.SLPmatVecMultiplyGalerkin(X,innerGeom),...
%      rhs2(:),prams.gmresTol,prams.maxIter);
%end
%% do SYMMLQ to find density function
%om.writeStars
%message = ['****    symmlq took ' num2str(toc,'%4.2e') ...
%    ' seconds     ****'];
%om.writeMessage(message,'%s\n');
%if (flag == 0)
%  message = ['****    symmlq required ' num2str(iter,'%3d'),...
%      ' iterations    ****'];
%  om.writeMessage(message,'%s\n');
%elseif (flag == 1)
%  message = '****    symmlq tolerance not achieved    ****';
%  om.writeMessage(message,'%s\n');
%  message = ['****    achieved tolerance is ' ...
%      num2str(relres,'%4.2e') '   ****'];
%  om.writeMessage(message,'%s\n');
%  message = ['****    symmlq took ' num2str(iter,'%3d'),...
%      ' iterations        ****'];
%  om.writeMessage(message,'%s\n');
%else 
%  message = 'SYMMLQ HAD A PROBLEM';
%  om.writeMessage(message,'%s\n');
%  flag
%end
%om.writeStars
%om.writeMessage(' ');
%
%eta1 = reshape(eta1,2*prams.Ninner,prams.nv);
%eta2 = reshape(eta2,2*prams.Ninner,prams.nv);
%eta2 = eta2./sqrt(2*pi*[sa;sa])*sqrt(innerGeom.N);
%eta3 = reshape(eta3,2*prams.Ninner,prams.nv);
%eta3 = eta3./sqrt(2*pi*[sa;sa])*sqrt(innerGeom.N);
%
%%norm(op.SLPmatVecMultiply(eta1(:),innerGeom,0) - rhs(:))/norm(rhs(:))
%%norm(op.SLPmatVecMultiply(eta2(:),innerGeom,0) - rhs(:))/norm(rhs(:))
%%norm(op.SLPmatVecMultiply(eta3(:),innerGeom,0) - rhs(:))/norm(rhs(:))
%%
%%
%relResVec = [];
%%relResVec = relresvec1;
%%relResVec(1:numel(relresvec2),2) = relresvec2;
%%relResVec(1:numel(relresvec3),3) = relresvec3;
%

S = [];
P = S;
%S = zeros(2*prams.Ninner*prams.nv);
%P = zeros(2*prams.Ninner*prams.nv);
%e = eye(2*prams.Ninner*prams.nv,1);
%for j = 1:2*prams.Ninner*prams.nv
%%  disp(2*prams.Ninner*prams.nv - j)
%  S(:,j) = op.SLPmatVecMultiply(e,innerGeom,0);
%  P(:,j) = op.matVecInvBD(e,innerGeom,[]);
%  e(j) = 0;
%  e(j+1) = 1;
%end
