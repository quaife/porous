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
rhs = [1/2*ones(prams.Ninner*prams.nv,1); ...
    ones(prams.Ninner*prams.nv,1)];
% right-hand side which corresponds to no-slip on the solid walls

op = poten(prams.Ninner,options.fmm);
% object for evaluating layer potentials

%tic
%[eta,flag,relres,iter,relresvec] = gmres(...
%    @(X) op.SLPmatVecMultiply(X,innerGeom),...
%    rhs,[],prams.gmresTol,prams.maxIter,...
%    @(X) op.matVecPreco(X,innerGeom));
%%[eta,flag,relres,iter,relresvec] = gmres(...
%%    @(X) op.SLPmatVecMultiply(X,innerGeom),...
%%    rhs,[],prams.gmresTol,prams.maxIter);
%% do unpreconditioned GMRES to find density function
%om.writeStars
%message = ['****    pGMRES took ' num2str(toc,'%4.2e') ...
%    ' seconds     ****'];
%om.writeMessage(message,'%s\n');
%if (flag == 0)
%  message = ['****    pGMRES required ' num2str(iter(2),'%3d'),...
%      ' iterations    ****'];
%  om.writeMessage(message,'%s\n');
%elseif (flag == 1)
%  message = '****    GMRES tolerance not achieved     ****';
%  om.writeMessage(message,'%s\n');
%  message = ['****    achieved tolerance is ' ...
%      num2str(relres,'%4.2e') '   ****'];
%  om.writeMessage(message,'%s\n');
%  message = ['****    pGMRES took ' num2str(iter(2),'%3d'),...
%      ' iterations        ****'];
%  om.writeMessage(message,'%s\n');
%else 
%  message = 'GMRES HAD A PROBLEM';
%  om.writeMessage(message,'%s\n');
%  relresvec
%end
%om.writeStars
%om.writeMessage(' ');


sa = innerGeom.sa;
rhs2 = rhs.*sqrt(2*pi*[sa(:);sa(:)])/sqrt(innerGeom.N);

[eta,flag,relres,iter,relresvec] = gmres(...
    @(X) op.SLPmatVecMultiply2(X,innerGeom),...
    rhs2,[],prams.gmresTol,prams.maxIter);

eta = eta*sqrt(innerGeom.N)./sqrt(2*pi*[sa(:);sa(:)]);



%tic
%[eta,flag,relres,iter,relresvec] = pcg(...
%    @(X) op.SLPmatVecMultiply(X,innerGeom),...
%    rhs,prams.gmresTol,prams.maxIter,...
%    @(X) op.matVecPreco(X,innerGeom));
%om.writeStars
%message = ['****    pcg took ' num2str(toc,'%4.2e') ...
%    ' seconds        ****'];
%om.writeMessage(message,'%s\n');
%if (flag == 0)
%  message = ['****    pGMRES required ' num2str(iter,'%3d'),...
%      ' iterations    ****'];
%  om.writeMessage(message,'%s\n');
%elseif (flag == 1)
%  message = '****    GMRES tolerance not achieved     ****';
%  om.writeMessage(message,'%s\n');
%  message = ['****    achieved tolerance is ' ...
%      num2str(relres,'%4.2e') '   ****'];
%  om.writeMessage(message,'%s\n');
%  message = ['****    pGMRES took ' num2str(iter,'%3d'),...
%      ' iterations        ****'];
%  om.writeMessage(message,'%s\n');
%else 
%  message = 'PCG HAD A PROBLEM';
%  om.writeMessage(message,'%s\n');
%  flag
%end
%om.writeStars
%om.writeMessage(' ');
%
%
S = zeros(2*prams.Ninner*prams.nv);

e = eye(2*prams.Ninner*prams.nv,1);
for k = 1:2*prams.Ninner*prams.nv;
  S(:,k) = op.SLPmatVecMultiply2(e,innerGeom);
  e(k) = 0;
  e(k+1) = 1;
end

%S = [];


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
