function S = SLPsolve(Xinner,options,prams)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% DESCRIPTION:
% INPUTS
% Xinner - parameterization of the inner boundaries
% prams - set of parameters
% options - set of options
% OUTPUTS:
% none 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

om = monitor(options,prams);
innerGeom = capsules(Xinner,'inner');
% create objects for the inner and outer boundaries.
% The outer boundary will need more points so need
% two classes.  Also, think that the double-layer potential
% may be more appropriate for the outer boundary as the
% circle may not be a very good preconditioner.  This way, we can use
% multigrid with Picard

%rhs = [zeros(prams.Ninner*prams.nv,1); ...
%    reshape(innerGeom.X(end/2+1:end,:),prams.Ninner*prams.nv,1)];
%rhs = [1/2*ones(prams.Ninner*prams.nv,1); ...
%    ones(prams.Ninner*prams.nv,1)];
rhs = [1/2*ones(prams.Ninner,prams.nv); ones(prams.Ninner,prams.nv)];
% right-hand side which corresponds to no-slip on the solid walls

op = poten(prams.Ninner,options.fmm);
% object for evaluating layer potentials

tic
[eta,flag,relres,iter,relresvec] = gmres(...
    @(X) op.SLPmatVecMultiply(X,innerGeom),...
    rhs(:),[],prams.gmresTol,prams.maxIter,...
    @(X) op.matVecPreco(X,innerGeom));
%[eta,flag,relres,iter,relresvec] = gmres(...
%    @(X) op.SLPmatVecMultiply(X,innerGeom),...
%    rhs,[],prams.gmresTol,prams.maxIter);
% do unpreconditioned GMRES to find density function
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
  relresvec
end
om.writeStars
om.writeMessage(' ');

sa = innerGeom.sa;
rhs2 = rhs.*sqrt(2*pi*[sa;sa])/sqrt(innerGeom.N);

tic
[eta2,flag,relres,iter,relresvec] = minres(...
    @(X) op.SLPmatVecMultiply2(X,innerGeom),...
    rhs2(:),prams.gmresTol,prams.maxIter,...
    @(X) op.matVecPreco(X,innerGeom));
%[eta2,flag,relres,iter,relresvec] = minres(...
%    @(X) op.SLPmatVecMultiply2(X,innerGeom),...
%    rhs2,prams.gmresTol,prams.maxIter);
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
[eta3,flag,relres,iter,relresvec] = symmlq(...
    @(X) op.SLPmatVecMultiply2(X,innerGeom),...
    rhs2(:),prams.gmresTol,prams.maxIter,...
    @(X) op.matVecPreco(X,innerGeom));
%[eta3,flag,relres,iter,relresvec] = minres(...
%    @(X) op.SLPmatVecMultiply2(X,innerGeom),...
%    rhs2,prams.gmresTol,prams.maxIter);
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


eta2 = reshape(eta2,2*prams.Ninner,prams.nv);
eta2 = eta2./sqrt(2*pi*[sa;sa])*sqrt(innerGeom.N);
eta2 = eta2(:);
norm(eta - eta2(:))
eta3 = reshape(eta3,2*prams.Ninner,prams.nv);
eta3 = eta3./sqrt(2*pi*[sa;sa])*sqrt(innerGeom.N);
eta3 = eta3(:);
norm(eta - eta3(:))





S = zeros(2*prams.Ninner*prams.nv);

%e = eye(2*prams.Ninner*prams.nv,1);
%for k = 1:2*prams.Ninner*prams.nv;
%  disp(2*prams.Ninner*prams.nv - k);
%  S(:,k) = op.SLPmatVecMultiply2(e,innerGeom);
%  e(k) = 0;
%  e(k+1) = 1;
%end



%sigmaInner = zeros(2*prams.Ninner,prams.nv);
%% initialize space for desnity function along inner boundaries
%for k = 1:prams.nv
%  istart = 2*(k-1)*prams.Ninner + 1;
%  iend = istart + 2*prams.Ninner - 1;
%  sigmaInner(:,k) = eta(istart:iend);
%end
%% unstack the density function at the inner boundaries
%
%
