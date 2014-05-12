function [S,P] = SLPsolve(Xinner,options,prams)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% DESCRIPTION:
% INPUTS
% Xinner - parameterization of the inner boundaries
% prams - set of parameters
% options - set of options
% OUTPUTS:
% none 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

<<<<<<< HEAD
preco = false
=======
preco = false;
>>>>>>> c23e11a9e6eb087926de1d86f5e48957af03f4d1
op = poten(prams.Ninner,options.fmm);
om = monitor(options,prams);
innerGeom = capsules(Xinner,'inner');
% create objects for the inner and outer boundaries.
% The outer boundary will need more points so need
% two classes.  Also, think that the double-layer potential
% may be more appropriate for the outer boundary as the
% circle may not be a very good preconditioner.  This way, we can use
% multigrid with Picard

%rhs = [1/2*ones(prams.Ninner,prams.nv); ones(prams.Ninner,prams.nv)];
theta = (0:prams.Ninner-1)'*2*pi/prams.Ninner;
etaTrue = zeros(2*prams.Ninner,prams.nv);
for k = 1:prams.nv
  etaTrue(:,k) = [exp(cos(theta));exp(sin(theta))];
end
rhs = op.SLPmatVecMultiply(etaTrue,innerGeom);
rhs = reshape(rhs,2*prams.Ninner,prams.nv);
% right-hand side which corresponds to no-slip on the solid walls

% object for evaluating layer potentials

tic
if preco
  fprintf('Using preconditioner\n')
  [eta1,flag,relres,iter,relresvec1] = gmres(...
      @(X) op.SLPmatVecMultiply(X,innerGeom),...
      rhs(:),[],prams.gmresTol,prams.maxIter,...
      @(X) op.matVecPreco(X,innerGeom));
else
  fprintf('Not using preconditioner\n')
  [eta1,flag,relres,iter,relresvec1] = gmres(...
      @(X) op.SLPmatVecMultiply(X,innerGeom),...
      rhs(:),[],prams.gmresTol,prams.maxIter);
end
% do GMRES to find density function
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




sa = innerGeom.sa;
rhs2 = rhs.*sqrt(2*pi*[sa;sa])/sqrt(innerGeom.N);

tic
if preco
  [eta2,flag,relres,iter,relresvec2] = minres(...
      @(X) op.SLPmatVecMultiply2(X,innerGeom),...
      rhs2(:),prams.gmresTol,prams.maxIter,...
      @(X) op.matVecPreco(X,innerGeom));
else
  [eta2,flag,relres,iter,relresvec2] = minres(...
      @(X) op.SLPmatVecMultiply2(X,innerGeom),...
      rhs2(:),prams.gmresTol,prams.maxIter);
end
% do MINRES to find density function
om.writeStars
message = ['****    minres took ' num2str(toc,'%4.2e') ...
    ' seconds     ****'];
om.writeMessage(message,'%s\n');
if (flag == 0)
  message = ['****    minres required ' num2str(iter,'%3d'),...
      ' iterations    ****'];
  om.writeMessage(message,'%s\n');
elseif (flag == 1)
  message = '****    minres tolerance not achieved    ****';
  om.writeMessage(message,'%s\n');
  message = ['****    achieved tolerance is ' ...
      num2str(relres,'%4.2e') '   ****'];
  om.writeMessage(message,'%s\n');
  message = ['****    minres took ' num2str(iter,'%3d'),...
      ' iterations        ****'];
  om.writeMessage(message,'%s\n');
else 
  message = 'MINRES HAD A PROBLEM';
  om.writeMessage(message,'%s\n');
  flag
end
om.writeStars
om.writeMessage(' ');


tic
if preco
  [eta3,flag,relres,iter,relresvec3] = symmlq(...
      @(X) op.SLPmatVecMultiply2(X,innerGeom),...
      rhs2(:),prams.gmresTol,prams.maxIter,...
      @(X) op.matVecPreco(X,innerGeom));
else
  [eta3,flag,relres,iter,relresvec3] = symmlq(...
      @(X) op.SLPmatVecMultiply2(X,innerGeom),...
      rhs2(:),prams.gmresTol,prams.maxIter);
end
% do SYMMLQ to find density function
om.writeStars
message = ['****    symmlq took ' num2str(toc,'%4.2e') ...
    ' seconds     ****'];
om.writeMessage(message,'%s\n');
if (flag == 0)
  message = ['****    symmlq required ' num2str(iter,'%3d'),...
      ' iterations    ****'];
  om.writeMessage(message,'%s\n');
elseif (flag == 1)
  message = '****    symmlq tolerance not achieved    ****';
  om.writeMessage(message,'%s\n');
  message = ['****    achieved tolerance is ' ...
      num2str(relres,'%4.2e') '   ****'];
  om.writeMessage(message,'%s\n');
  message = ['****    symmlq took ' num2str(iter,'%3d'),...
      ' iterations        ****'];
  om.writeMessage(message,'%s\n');
else 
  message = 'SYMMLQ HAD A PROBLEM';
  om.writeMessage(message,'%s\n');
  flag
end
om.writeStars
om.writeMessage(' ');

eta1 = reshape(eta1,2*prams.Ninner,prams.nv);
eta2 = reshape(eta2,2*prams.Ninner,prams.nv);
eta2 = eta2./sqrt(2*pi*[sa;sa])*sqrt(innerGeom.N);
eta3 = reshape(eta3,2*prams.Ninner,prams.nv);
eta3 = eta3./sqrt(2*pi*[sa;sa])*sqrt(innerGeom.N);

%norm(op.SLPmatVecMultiply(eta1(:),innerGeom) - rhs(:))/norm(rhs(:))
%norm(op.SLPmatVecMultiply(eta2(:),innerGeom) - rhs(:))/norm(rhs(:))
%norm(op.SLPmatVecMultiply(eta3(:),innerGeom) - rhs(:))/norm(rhs(:))
%
%
relResVec = [];
%relResVec = relresvec1;
%relResVec(1:numel(relresvec2),2) = relresvec2;
%relResVec(1:numel(relresvec3),3) = relresvec3;


S = [];
P = S;
S = zeros(2*prams.Ninner*prams.nv);
P = zeros(2*prams.Ninner*prams.nv);
e = eye(2*prams.Ninner*prams.nv,1);
for j = 1:2*prams.Ninner*prams.nv
%  disp(2*prams.Ninner*prams.nv - j)
  S(:,j) = op.SLPmatVecMultiply(e,innerGeom);
  P(:,j) = op.matVecPreco(e,innerGeom);
  e(j) = 0;
  e(j+1) = 1;
end
