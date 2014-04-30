classdef poten 
% this class defines single and double layers for various kernels 
% (stokes, laplace) on 2D periodic curves.  Also defines the 
% integrals required for pressure and stress.
% Defines the matricies that map a density function defined on the
% boundary of a curve to the layer potential evaluated on the 
% curve, and defines the operators that takes a density function 
% defined on the boundary of a curve and returns the layer 
% potential at arbitrary target points.
% This class also has the main routine that evaluates layer
% potentials using near-singular integration.
    
properties
  qw; 
  % quadrature weights for logarithmic singularity
  qp; 
  % quadrature points for logarithmic singularity (Alpert's rule)
  interpMat;  
  % interpolation matrix used for near-singular integration
  % This matrix replaces the need to use polyfit
  fmm
  % use the fmm or not
end % properties

methods 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function o = poten(N,fmm)
% o = poten(N,fmm): constructor; N is the number of points per curve.
% initialize class.

o.interpMat = o.lagrangeInterp;
% load in the interpolation matrix which is precomputed
% with 7 interpolation points
accuracyOrder = 8;
qw = o.quadratureS(N/2, accuracyOrder);
o.qp = qw(:,2:end);
o.qw = qw(:,1);
o.fmm = fmm;

end % poten: constructor
  

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function qw = quadratureS(o,m,q);
% qw = quadratureS(m,q) generates the quadrature rules for a function
% with m points and a logarithmic singularity at the origin.  q 
% controls the accuracy.  All rules are from Alpert 1999.

[v,u,a] = o.getWeights(q,2);
[x,w] = o.regularWeights(a);

n = m-2*a+1;
h = 1/m;

evalpots1 = v*h;
evalpots3 = 1-x*h;

ys = pi*[evalpots1;evalpots3];
wt = pi*h*[u;w];

wt = [wt; flipud(wt)]; ys = [ys; 2*pi-flipud(ys)];
oc = curve;
A = oc.sinterpS(2*m, ys); 
h = pi/m; 
yt = [a*h:h:(m - a)*h]';
%regular points away from singularity
wt = [wt; h*ones(2*length(yt),1)]/4/pi;
%quadrature weights at regular points away from singularity
lyt = 2*length(yt); 
B = sparse(lyt, 2*m); 
pos = 1+[(a:m-a)'; (m+a:2*m-a)'];

for k = 1:lyt
  B(k, pos(k)) = 1;
end
A = [sparse(A); B];
qw = [wt, A];



end % quadratureS

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [v,u,a] = getWeights(o,q,ker)
% [v,u,a] = getWeights(q,ker) loads quadrature rules for 
% different types of singularities.  All rules come from 
% Bradley Alpert's papers.  We are interested in nodesLogx.dat;

switch ker
case 0
  xp = load('nodesr.dat');
  lenth = [2;4;8];
  par = [2;4;7];
case 1
  xp = load('nodes_sqrtx.dat');
  lenth = [4;8;16];
  par = [3;5;10];
case 2
  xp = load('nodesLog.dat');
  lenth = [3;7;15];
  par = [2;5;10];
end

switch q
case 4
  v = xp(1:lenth(1), 1);
  u = xp(1:lenth(1), 2);
  a = par(1);
case 8
  v = xp(1+lenth(1):lenth(1)+lenth(2),1);
  u = xp(1+lenth(1):lenth(1)+lenth(2),2);
  a = par(2);
case 16
  v = xp(1+lenth(2)+lenth(1):sum(lenth),1);
  u = xp(1+lenth(2)+lenth(1):sum(lenth),2);
  a = par(3);
end


end % getWeights

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [x,w] = regularWeights(o,a)
% [x,w] = regularWeights(a) gets quadrature rules
% for regular functions (Alpert's quadrature rules).  These are
% needed for integrating away from singularities

par = [2;3;4;5;7;10];
for i = 1:length(par)
  if par(i)==a
    key = i;
  end
end

lenth = [2;3;4;6;8;12];
xp = load('nodesRegular.dat');
if key==1
  starting = 1;
else
  starting = sum(lenth(1:key-1))+1; 
end

x = xp( starting:starting+lenth(key)-1,1);
w = xp(starting:starting+lenth(key)-1,2);

end % regularWeights


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function invGf = matVecPreco(o,f,innerGeom,outerGeom)

invGf = f;

sigmah = zeros(2*innerGeom.N,1);
modes = (-innerGeom.N/2:innerGeom.N/2-1)';
for k = 1:innerGeom.nv
  rad = innerGeom.length(k)/2/pi;
  istart = (k-1)*2*innerGeom.N+1;
  iend = istart + 2*innerGeom.N - 1;
  sigma = f(istart:iend);
  sigmah(1:end/2) = fftshift(fft(sigma(1:end/2)));
  sigmah(end/2+1:end) = fftshift(fft(sigma(end/2+1:end)));
  for j = 1:innerGeom.N
    if abs(modes(j)) > 1
      % these terms are all diagonal in Fourier space
      const = rad/(4*abs(modes(j)));
      sigmah(j) = sigmah(j)/const;
      sigmah(j+innerGeom.N) = sigmah(j+innerGeom.N)/const;
    elseif abs(modes(j)) == 1
      % the -1 and 1 modes communicate, so save these to
      % compute preconditioner a bit later
      if modes(j) == -1
        fm1 = [sigmah(j);sigmah(j+innerGeom.N)];
      else
        fp1 = [sigmah(j);sigmah(j+innerGeom.N)];
      end
    else
      % 0 mode is diagonal in Fourier space
      const = -rad*log(rad)/2+rad/4;
      sigmah(j) = sigmah(j)/const;
      sigmah(j+innerGeom.N) = sigmah(j+innerGeom.N)/const;
    end

  end

  A = 1/2/rad*[7 0 1 0; 0 7 0 -1; 1 0 7 0; 0 -1 0 7] + ...
     1i/2/rad*[0 1 0 1; -1 0 1 0; 0 -1 0 -1; -1 0 1 0];
% inverse of the matrix does the following four mappings for functions
% that only contain 1 and -1 modes
%    [sin(theta);0]           -> r/8*[3*sin(theta);-cos(theta)] 
%    [0;cos(theta)]           -> r/8*[-sin(theta);3*cos(theta)]
%    [cos(theta);-sin(theta)] -> r/4*[cos(theta);-sin(theta)]
%    [cos(theta);sin(theta)]  -> r/4*[cos(theta);sin(theta)]
% The image of [cos(theta);sin(theta)] under the SLP is in fact 0, but
% we do the above mapping to make the preconditioner invertible and to
% keep this mode consistent with the [cos(theta);-sin(theta)] mapping

  g = A*[fm1;fp1]; 
  j = innerGeom.N/2;
  sigmah(j) = g(1);
  sigmah(j+innerGeom.N) = g(2);
  j = innerGeom.N/2+2;
  sigmah(j) = g(3);
  sigmah(j+innerGeom.N) = g(4);


  invGf(istart:iend) = ...
      [ifft(ifftshift(sigmah(1:end/2))); ...
       ifft(ifftshift(sigmah(end/2+1:end)))];

end % k = exclusions



invGf = real(invGf);
%invGf = f;



end % matVecPreco

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function Gf = matVecMultiply(o,f,innerGeom,outerGeom)

Ninner = innerGeom.N;
nv = innerGeom.nv;
Nouter = outerGeom.N;
Gfinner = zeros(2*Ninner,nv);
Gfouter = zeros(2*Nouter,1);
innerEta = zeros(2*Ninner,nv);
outerEta = zeros(2*Nouter,1);

for k = 1:nv
  istart = (k-1)*2*Ninner + 1;
  iend = istart + 2*Ninner - 1;
  innerEta(:,k) = f(istart:iend);
end
istart = nv*2*Ninner + 1;
iend = istart + 2*Nouter - 1;
outerEta = f(istart:iend);
% unstack f so that it is one x-coordinate and one y-coordinate per
% column


Gfinner = Gfinner + o.exactStokesSLdiag(innerGeom,innerEta);
% diagonal term from exclusions

Gfouter = Gfouter - 0.5*outerEta;
% jump term
Gfouter = Gfouter + o.exactStokesDLdiag(outerGeom,outerEta);
% double-layer potential
Gfouter = Gfouter + o.exactStokesN0diag(outerGeom,outerEta);
% rank one modification to remove null space

if ~o.fmm
  stokesSLP = exactStokesSL(o,innerGeom,innerEta);
  [~,stokesSLPtar] = ...
      exactStokesSL(o,innerGeom,innerEta,outerGeom.X,(1:nv));
else
  stokesSLP = exactStokesSLfmm(o,innerGeom,innerEta);
  [~,stokesSLPtar] = ...
      exactStokesSLfmm(o,innerGeom,innerEta,outerGeom.X,(1:nv));
end

Gfinner = Gfinner + stokesSLP;
% add in contribution from all other exclusions
Gfouter = Gfouter + stokesSLPtar;
% add in contribution from all exclusions to outer boundary

[~,stokesDLPtar] = ...
    exactStokesDL(o,outerGeom,outerEta,innerGeom.X,1);
Gfinner = Gfinner + stokesDLPtar;
% add in contribution from the outer boundary to all the exclusions
    

Gf = [Gfinner(:);Gfouter(:)];


end % matVecMultiply



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function stokesSLP = exactStokesSLdiagOld(o,geom,f)

X = geom.X;
N = size(X,1)/2;
oc = curve;
[x,y] = oc.getXY(X);
sa = geom.sa;
[fx,fy] = oc.getXY(f.*[sa;sa]);
stokesSLP = zeros(2*N,1);

for j = 1:N
  ind = 1 + mod(j-1 + (0:N-1),N);
  xSou = qp*x(ind);
  ySou = qp*y(ind);
  fxSou = qp*fx(ind);
  fySou = qp*fy(ind);

  rho = (x(j) - xSou).^2 + (y(j) - ySou).^2;
  br = 1./rho;
  logpart = -1/2*qw.*log(rho);
  stokesSLP(j) = sum((logpart.*fxSou)'*qp);
  stokesSLP(j+N) = sum((logpart.*fySou)'*qp);
  % stokes single-layer potential due to the log componenet

  rdotf = qw.*((x(j) - xSou).*fxSou + (y(j) - ySou).*fySou).*br;
  stokesSLP(j) = stokesSLP(j) + sum((rdotf .* (x(j) - xSou))'*qp);
  stokesSLP(j+N) = stokesSLP(j+N) + sum((rdotf .* (y(j) - ySou))'*qp);
end

end % exactStokesSLdiag




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function stokesSLP = exactStokesSLdiag(o,geom,f)

X = geom.X;
N = size(X,1)/2;
nv = size(X,2);
oc = curve;
[x,y] = oc.getXY(X);
sa = geom.sa;
[fx,fy] = oc.getXY(f.*[sa;sa]);
stokesSLP = zeros(2*N,nv);
qw = o.qw;
qp = o.qp;

for k = 1:nv
  for j = 1:N
    ind = 1 + mod(j-1 + (0:N-1),N);
    xSou = qp*x(ind,k);
    ySou = qp*y(ind,k);
    fxSou = qp*fx(ind,k);
    fySou = qp*fy(ind,k);

    rho = (x(j,k) - xSou).^2 + (y(j,k) - ySou).^2;
    br = 1./rho;
    logpart = -1/2*qw.*log(rho);
    stokesSLP(j,k) = sum((logpart.*fxSou)'*qp);
    stokesSLP(j+N,k) = sum((logpart.*fySou)'*qp);
    % stokes single-layer potential due to the log componenet

    rdotf = qw.*((x(j,k) - xSou).*fxSou + (y(j,k) - ySou).*fySou).*br;
    stokesSLP(j,k) = stokesSLP(j,k) + sum((rdotf .* (x(j,k) - xSou))'*qp);
    stokesSLP(j+N,k) = stokesSLP(j+N,k) + sum((rdotf .* (y(j,k) - ySou))'*qp);
  end
end

end % exactStokesSLdiag



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function G = stokesSLmatrix(o,vesicle)
% G = stokesSLmatrix(vesicle) generate single layer potential for
% Stokes vesicle is a data structure defined as in the curve class
% G is (2N,2N,nv) array where N is the number of points per curve and 
% nv is the number of curves in X 

X = vesicle.X; % position of vesicle
oc = curve;
[x,y] = oc.getXY(X);
% seperate x and y coordinates of vesicle
sa = vesicle.sa; % Jacobian
N = size(X,1)/2;

G = zeros(2*N,2*N,vesicle.nv);
% initalize single-layer potential to zero
for k=1:vesicle.nv  % Loop over curves
  for j=1:N % Loop over targets
    ind = 1 + mod(j-1 + (0:N-1),N);
    xin = o.qp*x(ind,k);
    yin = o.qp*y(ind,k);

    rho = (xin-x(j,k)).^2 + (yin-y(j,k)).^2;
    br = 1./rho;

    logpart = -1/2 * o.qw .* log(rho);
    
    G(j, ind,   k) = (logpart+(o.qw.*...
        (( x(j,k) - xin).^2.*br)))'*o.qp;
    G(j, N+ind, k) = (o.qw.* ...
        (x(j,k)-xin).*(y(j,k)-yin).*br)'*o.qp;
    G(j+N, 1:N ,k) = G(j,N+1:2*N,k);
    G(j+N,N+ind,k) = (logpart+(o.qw.*...
        ( ( y(j,k) - yin).^2.*br)))'*o.qp;
  end % j
  sak = repmat(sa(:,k)',2*N,2);
  
  G(:,:,k) = G(:,:,k).*sak;
  % multiply both components of G by the Jacobian
end % k


end % stokesSLmatrix



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function stokesDLP = exactStokesDLdiag(o,geom,f)

oc = curve;
[x,y] = oc.getXY(geom.X);
[nx,ny] = oc.getXY(geom.normal);
[tx,ty] = oc.getXY(geom.xt);
[fx,fy] = oc.getXY(f);
N = size(x,1);
sa = geom.sa;
fDOTt = fx.*tx + fy.*ty;

stokesDLP = zeros(2*N,1);


for j = 1:N
  rho2 = (x(j)-x).^2 + (y(j)-y).^2;

  coeff = ((x(j) - x).*nx + (y(j) - y).*ny).*...
          ((x(j) - x).*fx + (y(j) - y).*fy).*...   
          sa./rho2.^2/pi*2*pi/N;
  coeff(j) = 0;
  % Set diagonal term to one to avoid dividing by zero

  stokesDLP(j) = sum(coeff.*(x(j) - x));
  stokesDLP(j+N) = sum(coeff.*(y(j) - y));
  stokesDLP(j) = stokesDLP(j) - ...
      1/2/pi*geom.cur(j)*fDOTt(j)*tx(j)*sa(j)*2*pi/N;
  stokesDLP(j+N) = stokesDLP(j+N) - ...
      1/2/pi*geom.cur(j)*fDOTt(j)*ty(j)*sa(j)*2*pi/N;
end

end % exactStokesDLdiag


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function D = stokesDLmatrix(o,vesicle)
% D = stokesDLmatrix(vesicle), generate double-layer potential for 
% Stokes vesicle is a data structure defined as in the capsules class
% D is (2N,2N,nv) array where N is the number of points per curve and 
% nv is the number of curves in X 


N = vesicle.N; % Number of points per vesicle
nv = vesicle.nv; % Number of vesicles
X = vesicle.X; % Vesicle positions
oc = curve;
[x,y] = oc.getXY(X);

normal = vesicle.normal; % Normal vector
tang = vesicle.xt; % Tangent vector
sa = vesicle.sa; % Jacobian

D = zeros(2*N,2*N,nv);
for k=1:nv  % Loop over curves
  for j=1:N % Loop over targets
    rho2 = (x(j,k)-x(:,k)).^2 + (y(j,k)-y(:,k)).^2;
    rho2(j) = 1;
    % Set diagonal term to one to avoid dividing by zero

    coeff = ((x(j,k) - x(:,k)).*normal(1:N,k) + ...
        (y(j,k) - y(:,k)).*normal(N+1:2*N,k)).*...
        sa(:,k)./rho2.^2/pi;
    % part of kernel of the double-layer potential

    D(j,:,k) = 2*pi/N*[coeff.*(x(j,k) - x(:,k)).^2; ...
      coeff.*(x(j,k) - x(:,k)).*(y(j,k) - y(:,k))]';
    D(j+N,:,k) = 2*pi/N*[coeff.*(y(j,k) - y(:,k)).*(x(j,k)-x(:,k)); ...
      coeff.*(y(j,k) - y(:,k)).^2]';
    % Build double-layer potential matrix D

    rc = [j j+N];
    D(rc,rc,k) = -2*pi/N*sa(j,k)*vesicle.cur(j,k)/2/pi*...
      [tang(j,k)^2 tang(j,k)*tang(j+N,k);...
      tang(j+N,k)*tang(j,k) tang(j+N,k)^2];
    % Diagonal term requires the limiting value.  Otherwise, above formula
    % would divide by zero
  end % j
end % k


end % stokesDLmatrix



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function stokesN0 = exactStokesN0diag(o,geom,f)

oc = curve;
[nx,ny] = oc.getXY(geom.normal);
[fx,fy] = oc.getXY(f);
N = size(nx,1);
sa = geom.sa;
INTnDOTf = sum((nx.*fx + ny.*fy).*sa)*2*pi/N;

stokesN0 = zeros(2*N,1);

for j = 1:N
  stokesN0(j) = INTnDOTf*nx(j);
  stokesN0(j+N) = INTnDOTf*ny(j);
end

end % exactStokesN0diag

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function N0 = stokesN0matrix(o,vesicle)
% N0 = stokesN0matrix(vesicle) generates the the integral 
% operator with kernel normal(x) \otimes normal(y) which removes
% the rank one defficiency of the double-layer potential.  Need
% this operator for solid walls

N = vesicle.N; % Number of points per vesicle
nv = vesicle.nv; % Number of vesicles
X = vesicle.X; % Vesicle positions
oc = curve;
[x,y] = oc.getXY(X);

normal = vesicle.normal; % Normal vector
sa = vesicle.sa; % Jacobian


N0 = zeros(2*N,2*N,nv);
for k = 1:1 % Loop over curves
  % Only want to form N0 for the outer boundary.
  % Rotlets and Stokeslets take care of null space for inner 
  % boundaries
  for j = 1:2*N % Loop over targets
    N0(j,:,k) = normal(j,k)*normal(:,k).*...
        [sa(:,k);sa(:,k)]*2*pi/N;
    % compute the columns of the modification to the double-layer
    % potential.
  end
end
% Use N0 if solving (-1/2 + DLP)\eta = f where f has no flux through the
% boundary.  By solving (-1/2 + DLP + N0)\eta = f, we guarantee that \eta
% also has no flux through the boundary.  This is not required, but it means
% we're enforcing one addition condition on eta which removes the rank
% one kernel.  DLP is the double-layer potential for stokes equation

end % stokesN0matrix

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function G = laplaceSLmatrix(o,vesicle)
% G = laplaceSLmatrix(vesicle), generate single layer potential for
% laplace.  vesicle is a data structure defined as in the capsules class
% G is (N,N,nv) array where N is the number of points per curve and 
% nv is the number of curves in X

N = vesicle.N; % Number of points per vesicle
nv = vesicle.nv; % Number of vesicles
X = vesicle.X; % Vesicle positions
oc = curve;
[x,y] = oc.getXY(X);

sa = vesicle.sa; % Jacobian

G = zeros(o.N,o.N,vesicle.nv);
% initalize single-layer potential to zero
for k=1:vesicle.nv  % Loop over curves
  for j=1:o.N % Loop over targets
    ind = 1 + mod(j-1 + (0:o.N-1),o.N);
    xin = o.qp*x(ind,k);
    yin = o.qp*y(ind,k);

    rho = (xin-x(j,k)).^2 + (yin-y(j,k)).^2;
    br = 1./rho;

    logpart = -1/2 * o.qw .* log(rho);
    G(j,ind,k) = logpart'*o.qp;
    
    G(j, ind,   k) = (logpart+(o.qw.*( ( x(j,k) - xin).^2 .* br)))'*o.qp;
    G(j, o.N+ind, k) = (o.qw.*(x(j,k)-xin).*(y(j,k)-yin) .* br)'*o.qp;
    G(j+o.N, 1:o.N ,k) = G(j,o.N+1:2*o.N,k);
    G(j+o.N,o.N+ind,k) = (logpart+(o.qw.*( ( y(j,k) - yin).^2 .* br)))'*o.qp;
  end % j
  sak = repmat(sa(:,k)',o.N,1);
  
  G(:,:,k) = G(:,:,k).*sak;
  % multiply G by the Jacobian
end % k


end % stokesSLmatrix



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function D = laplaceDLmatrix(o,vesicle)
% D = laplaceDLmatrix(vesicle), generate double layer potential for
% laplace.  vesicle is a data structure defined as in the capsules class
% D is (N,N,nv) array where N is the number of points per curve and 
% nv is the number of curves in X.  Need this for collision detection

N = vesicle.N; % Number of points per vesicle
nv = vesicle.nv; % Number of vesicles
X = vesicle.X; % Vesicle positions
oc = curve;
[x,y] = oc.getXY(X);

normal = vesicle.normal; % Normal vector
tang = vesicle.xt; % Tangent vector
sa = vesicle.sa; % Jacobian
cur = vesicle.cur; % curvature

D = zeros(N,N,nv);
% initialize double-layer potential to zero
for k=1:nv  % Loop over curves
  for j=1:N % Loop over targets
    rho2 = (x(:,k)-x(j,k)).^2 + (y(:,k)-y(j,k)).^2;
    rho2(j) = 1;
    % Set diagonal term to one to avoid dividing by zero

    coeff = ((x(:,k) - x(j,k)).*normal(1:N,k) + ...
        (y(:,k) - y(j,k)).*normal(N+1:2*N,k)).*...
        sa(:,k)./rho2/2/pi;
    % kernel of double-layer potential for Laplace
    coeff(j) = 0.5*cur(j,k)/2/pi*sa(j,k);
    % diagonal term is proportional to the curvature
    D(j,:,k) = 2*pi/N*coeff;
    % multiply by d\theta term

  end % j
end % k


end % laplaceDLmatrix


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function vel = velEval(o,t,Xtar,innerGeom,outerGeom,sigmaInner,sigmaOuter)

targetPnts = capsules(Xtar,'targets');

[~,NearInner] = innerGeom.getZone(targetPnts,2);
[~,NearOuter] = outerGeom.getZone(targetPnts,2);

if ~o.fmm
  vel1 = o.nearSingInt(innerGeom,sigmaInner,@o.exactStokesSLdiag,...
      NearInner,@o.exactStokesSL,targetPnts,0,'inner');
  vel2 = o.nearSingInt(outerGeom,sigmaOuter,@o.exactStokesDLdiag,...
      NearOuter,@o.exactStokesDL,targetPnts,0,'outer');
else
  vel1 = o.nearSingInt(innerGeom,sigmaInner,@o.exactStokesSLdiag,...
      NearInner,@o.exactStokesSLfmm,targetPnts,0,'inner');
  vel2 = o.nearSingInt(outerGeom,sigmaOuter,@o.exactStokesDLdiag,...
      NearOuter,@o.exactStokesDL,targetPnts,0,'outer');
end


vel = vel1 + vel2;

%min(min(Xtar(end/2:end,:)))
s = find(Xtar(end/2+1:end,:) < -5);
%fprintf('Current time is %4.2e\n',t);
%[t numel(s)]
vel(s) = 0;
vel(s+size(Xtar,1)/2) = 0;
%if numel(s) > 0
%  disp(s)
%  disp([Xtar vel])
%  pause
%end




end % velEval



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function LP = nearSingInt(o,souPts,f,diagLP,...
    NearStruct,kernel,tarPts,tEqualS,side,trash)
% LP = nearSingInt(souPts,f,diagLP,...
% zone,dist,nearest,icp,argnear,tarPts,tequalS) computes a 
% layer potential due to f at all points in tarPts.X.  If
% tEqualS==true, then the tarPts == souPts and the
% self-interaction is skipped.
% diagLP is the diagonal of the potential needed to compute
% the layer potential of each source curve indepenedent of all others
% if using near-singular integration
% zone,dist,nearest,icp,argnear are required by near-singular
% integration (they keep everything sorted and precomputed)
% Everything is in the 2*N x nv format
% Can pass a final argument if desired so that plots of the 
% near-singular integration algorithm are displayed

dist = NearStruct.dist;
zone = NearStruct.zone;
nearest = NearStruct.nearest;
icp = NearStruct.icp;
argnear = NearStruct.argnear;

Xsou = souPts.X; % source positions
Nsou = size(Xsou,1)/2; % number of source points
nvSou = size(Xsou,2); % number of source curves
Xtar = tarPts.X; % target positions
Ntar = size(Xtar,1)/2; % number of target points
nvTar = size(Xtar,2); % number of target curves

h = souPts.length/Nsou; % arclength term

%Nup = 2^ceil(3/2*log2(Nsou));
Nup = Nsou*2^ceil(1/2*log2(Nsou));
% upsample to N^(3/2).  
% only want to add on powers of 2 so that ffts are simple
% Nup at least has to be a multiple of N

vself = diagLP(souPts,f);
if strcmp(side,'outer')
  vself = vself - 0.5*f;
end
% Compute velocity due to each curve independent of others.
% This is needed when using near-singular integration since
% we require a value for the layer-potential on the curve of 
% sources 

Xup = [interpft(Xsou(1:Nsou,:),Nup);interpft(Xsou(Nsou+1:2*Nsou,:),Nup)];
fup = [interpft(f(1:Nsou,:),Nup);interpft(f(Nsou+1:2*Nsou,:),Nup)];
% upsample positions, traction jump

souUpPts = capsules(Xup,side);
% Build an object with the upsampled curve

interpOrder = size(o.interpMat,1);
% lagrange interpolation order
p = ceil((interpOrder+1)/2);
% want half of the lagrange interpolation points to the left of the
% closest point and the other half to the right
vel = zeros(2*Ntar,nvTar,nvSou);
% allocate space for storing velocity at intermediate points
% needed by near-singular integration

for k1 = 1:nvSou
  if tEqualS % sources == targets
    K = [(1:k1-1) (k1+1:nvTar)];
    % skip diagonal curve
  else % sources ~= targets
    K = (1:nvTar);
    % consider all curves
  end
  for k2 = K
    J = find(zone{k1}(:,k2) == 1);
    % set of points on curve k2 close to curve k1
    indcp = icp{k1}(J,k2);
    % closest point on curve k1 to each point on curve k2 
    % that is close to curve k1
    for j = 1:numel(J)
      pn = mod((indcp(j)-p+1:indcp(j)-p+interpOrder)' - 1,Nsou) + 1;
      % index of points to the left and right of the closest point
      v = filter(1,[1 -full(argnear{k1}(J(j),k2))],...
        o.interpMat*vself(pn,k1));
      vel(J(j),k2,k1) = v(end);
      % x-component of the velocity at the closest point
      v = filter(1,[1 -full(argnear{k1}(J(j),k2))],...
        o.interpMat*vself(pn+Nsou,k1));
      vel(J(j)+Ntar,k2,k1) = v(end);
      % y-component of the velocity at the closest point
    end
  end
end
% compute values of velocity at required intermediate points
% using local interpolant

if tEqualS % sources == targets
  farField = kernel(souUpPts,fup);
  % evaluate layer potential at all targets except ignore the
  % diagonal term
else % sources ~= targets
  [~,farField] = kernel(souUpPts,fup,Xtar,1:nvSou);
  % evaluate layer potential due to all curves at all points
  % in Xtar;
end
% Use upsampled trapezoid rule to compute layer potential

nearField = zeros(2*Ntar,nvTar);
% Initialize potential at near points to zero

beta = 1.1;
% small buffer to make sure Lagrange interpolation points are
% not in the near zone
for k1 = 1:nvSou
  if tEqualS % sources == targets
    K = [(1:k1-1) (k1+1:nvTar)];
    % skip diagonal curve
  else % sources ~= targets
    K = (1:nvTar);
    % consider all curves
  end
  for k2 = K
    J = find(zone{k1}(:,k2) == 1);
    if (numel(J) ~= 0)
      [~,potTar] = kernel(souUpPts,fup,...
          [Xtar(J,k2);Xtar(J+Ntar,k2)],k1);
      % Need to subtract off contribution due to curve k1 since
      % its layer potential will be evaulted using Lagrange
      % interpolant of nearby points
      nearField(J,k2) =  - potTar(1:numel(J));
      nearField(J+Ntar,k2) =  - potTar(numel(J)+1:end);

      XLag = zeros(2*numel(J),interpOrder - 1);
      for i = 1:numel(J)
        nx = (Xtar(J(i),k2) - nearest{k1}(J(i),k2))/...
            dist{k1}(J(i),k2);
        ny = (Xtar(J(i)+Ntar,k2) - nearest{k1}(J(i)+Ntar,k2))/...
            dist{k1}(J(i),k2);

        XLag(i,:) = nearest{k1}(J(i),k2) + ...
            beta*h(k2)*nx*(1:interpOrder-1);
        XLag(i+numel(J),:) = nearest{k1}(J(i)+Ntar,k2) + ...
            beta*h(k2)*ny*(1:interpOrder-1);
        % Lagrange interpolation points coming off of curve k1
        % All points are behind Xtar(J(i),k2) and are sufficiently
        % far from curve k1 so that the Nup-trapezoid rule gives
        % sufficient accuracy
      end
      [~,lagrangePts] = kernel(souUpPts,fup,XLag,k1);
      % evaluate velocity at the lagrange interpolation points
      x = XLag(1,:);
      y = XLag(2,:);

      for i = 1:numel(J)
        Px = o.interpMat*[vel(J(i),k2,k1) ...
            lagrangePts(i,:)]';
        Py = o.interpMat*[vel(J(i)+Ntar,k2,k1) ...
            lagrangePts(i+numel(J),:)]';
        % Build polynomial interpolant along the one-dimensional
        % points coming out of the curve 
        dscaled = full(dist{k1}(J(i),k2)/...
            (beta*h(k2)*(interpOrder-1)));
        % Point where interpolant needs to be evaluated

        v = filter(1,[1 -dscaled],Px);
        v1 = v;
        nearField(J(i),k2) = nearField(J(i),k2) + ...
            v(end);
        v = filter(1,[1 -dscaled],Py);
        v2 = v;
        nearField(J(i)+Ntar,k2) = nearField(J(i)+Ntar,k2) + ...
            v(end);
        % Assign higher-order results coming from Lagrange 
        % integration to velocity at near point.  Filter is faster
        % than polyval

        if (nargin == 10)
          figure(2); clf; hold on;
          plot(Xsou(1:Nsou,:),Xsou(Nsou+1:end,:),'r.')
%          plot(Xtar(1:Ntar,:),Xtar(Ntar+1:end,:),'k.')
          plot(Xtar(J,k2),Xtar(Ntar+J,k2),'b.')
          plot(XLag(1:numel(J),:),XLag(numel(J)+1:end,:),'kx')
          plot(XLag(i,:),XLag(numel(J)+i,:),'gx')
%          axis([3.5 4.5 23.5 24.5])
%          axis equal

          figure(1); clf; hold on
          plot((0:interpOrder-1)*beta*h(k2),...
              [vel(J(i),k2,k1) lagrangePts(i,:)],'g-o')
          plot((0:interpOrder-1)*beta*h(k2),...
              [vel(J(i)+Ntar,k2,k1) lagrangePts(i+numel(J),:)],'r--o')
          pause(0.01)
        end
        % DEBUG: PASS IN A DUMMY VARIABLE INTO THIS ROUTINE AND THEN
        % YOU CAN SEE THE INTERPOLATION POINTS AND CHECK THE SMOOTHNESS
        % OF THE INTERPOLANT

      end % i
    end % numel(J) ~= 0
    % Evaluate layer potential at Lagrange interpolation
    % points if there are any
  end % k2
end % k1

if tEqualS % souPts == tarPts
  LP = farField(1:Nup/Ntar:end,:) + nearField;
else % souPts ~= tarPts
  LP = farField + nearField;
end
% Add kernel due to far points and near points.  Far points were
% upsampled if source==target so need to truncate here.  We are 
% only using Ntar target points.  Note that it is only the sources 
% that were upsampled


end % nearSingInt



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [stokesSLP,stokesSLPtar] = ...
    exactStokesSL(o,vesicle,f,Xtar,K1)
% [stokesSLP,stokesSLPtar] = exactStokesSL(vesicle,f,Xtar,K1)
% computes the single-layer potential due to f around all vesicles 
% except itself.  Also can pass a set of target points Xtar and a 
% collection of vesicles K1 and the single-layer potential due to
% vesicles in K1 will be evaluated at Xtar.
% Everything but Xtar is in the 2*N x nv format
% Xtar is in the 2*Ntar x ncol format

X = vesicle.X; % Vesicle positions
sa = vesicle.sa; % Arclength term

if nargin == 5
  Ntar = size(Xtar,1)/2;
  ncol = size(Xtar,2);
  stokesSLPtar = zeros(2*Ntar,ncol);
else
  K1 = [];
  stokesSLPtar = [];
  ncol = 0;
  % if nargin ~= 5, the user does not need  the velocity at arbitrary
  % points
end

den = f.*[sa;sa]*2*pi/vesicle.N;
% multiply by arclength term

for k2 = 1:ncol % loop over columns of target points
  for j = 1:Ntar % loop over rows of target points
    dis2 = (Xtar(j,k2) - X(1:vesicle.N,K1)).^2 + ... 
        (Xtar(j+Ntar,k2) - X(vesicle.N+1:2*vesicle.N,K1)).^2;
    diffxy = [Xtar(j,k2) - X(1:vesicle.N,K1) ; ...
        Xtar(j+Ntar,k2) - X(vesicle.N+1:2*vesicle.N,K1)];
    % distance squared and difference of source and target location

    coeff = log(dis2)/2;
    % first part of single-layer potential for Stokes
    
    val = coeff.*den(1:vesicle.N,K1);
    stokesSLPtar(j,k2) = -sum(val(:));
    val = coeff.*den(vesicle.N+1:2*vesicle.N,K1);
    stokesSLPtar(j+Ntar,k2) = -sum(val(:));
    % log term in stokes single-layer potential

    coeff = (diffxy(1:vesicle.N,:).*den(1:vesicle.N,K1) + ...
        diffxy(vesicle.N+1:2*vesicle.N,:).*...
        den(vesicle.N+1:2*vesicle.N,K1))./dis2;
    % second part of single-layer potential for Stokes

    val = coeff.*diffxy(1:vesicle.N,:);
    stokesSLPtar(j,k2) = stokesSLPtar(j,k2) + sum(val(:));
    val = coeff.*diffxy(vesicle.N+1:2*vesicle.N,:);
    stokesSLPtar(j+Ntar,k2) = stokesSLPtar(j+Ntar,k2) + sum(val(:));
    % r \otimes r term of the stokes single-layer potential
  end % j

end % k2
% Evaluate single-layer potential at arbitrary target points
stokesSLPtar = 1/(4*pi)*stokesSLPtar;
% 1/4/pi is the coefficient in front of the single-layer potential

stokesSLP = zeros(2*vesicle.N,vesicle.nv); % Initialize to zero
if nargin == 3
  for k1 = 1:vesicle.nv % vesicle of targets
    K = [(1:k1-1) (k1+1:vesicle.nv)];
    % Loop over all vesicles except k1
    for j=1:vesicle.N
      dis2 = (X(j,k1) - X(1:vesicle.N,K)).^2 + ...
          (X(j+vesicle.N,k1) - X(vesicle.N+1:2*vesicle.N,K)).^2;
      diffxy = [X(j,k1) - X(1:vesicle.N,K) ; ...
          X(j+vesicle.N,k1) - X(vesicle.N+1:2*vesicle.N,K)];
      % distance squared and difference of source and target location

      coeff = log(dis2)/2;
      % first part of single-layer potential for Stokes

      val = coeff.*den(1:vesicle.N,K);
      stokesSLP(j,k1) = -sum(val(:));
      val = coeff.*den(vesicle.N+1:2*vesicle.N,K);
      stokesSLP(j+vesicle.N,k1) = -sum(val(:));
      % logarithm terms in the single-layer potential

      coeff = (diffxy(1:vesicle.N,:).*den(1:vesicle.N,K) + ...
          diffxy(vesicle.N+1:2*vesicle.N,:).*...
          den(vesicle.N+1:2*vesicle.N,K))./dis2;
      % second part of single-layer potential for Stokes

      val = coeff.*diffxy(1:vesicle.N,:);
      stokesSLP(j,k1) = stokesSLP(j,k1) + sum(val(:));
      val = coeff.*diffxy(vesicle.N+1:2*vesicle.N,:);
      stokesSLP(j+vesicle.N,k1) = stokesSLP(j+vesicle.N,k1) + ...
          sum(val(:));
      % r \otimes r term of the single-layer potential

    end % j
  end % k1
  % Evaluate single-layer potential at vesicles but oneself
  stokesSLP = 1/(4*pi)*stokesSLP;
  % 1/4/pi is the coefficient in front of the single-layer potential
end % nargin == 3

end % exactStokesSL

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [stokesDLP,stokesDLPtar] = ...
    exactStokesDL(o,vesicle,f,Xtar,K1)
% [stokesDLP,stokesDLPtar] = exactStokesDL(vesicle,f,Xtar,K1)
% computes the double-layer potential due to f around all vesicles 
% except itself.  Also can pass a set of target points Xtar and a 
% collection of vesicles K1 and the double-layer potential due to
% vesicles in K1 will be evaluated at Xtar.
% Everything but Xtar is in the 2*N x nv format
% Xtar is in the 2*Ntar x ncol format

nv = vesicle.nv; % number of vesicles
X = vesicle.X; % Vesicle positions
normal = vesicle.normal; % Outward normal
sa = vesicle.sa; % Jacobian

if nargin == 5
  Ntar = size(Xtar,1)/2;
  ncol = size(Xtar,2);
  stokesDLPtar = zeros(2*Ntar,ncol);
else
  K1 = [];
  stokesDLPtar = [];
  ncol = 0;
  Ntar = 0;
  % if nargin ~= 5, the user does not need the velocity at arbitrary
  % points
end

den = f.*[sa;sa]*2*pi/vesicle.N;

for k2 = 1:ncol % loop over columns of target points
  for j = 1:Ntar % loop over rows of target points
    diffxy = [Xtar(j,k2) - X(1:vesicle.N,K1) ; ...
        Xtar(j+Ntar,k2) - X(vesicle.N+1:2*vesicle.N,K1)];
    dis2 = diffxy(1:vesicle.N,:).^2 + ...
        diffxy(vesicle.N+1:2*vesicle.N,:).^2;
    % difference of source and target location and distance squared

    coeff = (diffxy(1:vesicle.N,:).*normal(1:vesicle.N,K1) + ...
      diffxy(vesicle.N+1:2*vesicle.N,:).*...
      normal(vesicle.N+1:2*vesicle.N,K1))./dis2.^2.* ...
      (diffxy(1:vesicle.N,:).*den(1:vesicle.N,K1) + ...
      diffxy(vesicle.N+1:2*vesicle.N,:).*...
      den(vesicle.N+1:2*vesicle.N,K1));
    % \frac{(r \dot n)(r \dot density)}{\rho^{4}} term

    stokesDLPtar(j,k2) = stokesDLPtar(j,k2) + ...
        sum(sum(coeff.*diffxy(1:vesicle.N,:)));
    stokesDLPtar(j+Ntar,k2) = stokesDLPtar(j+Ntar,k2) + ...
        sum(sum(coeff.*diffxy(vesicle.N+1:2*vesicle.N,:)));
    % r \otimes r term of the single-layer potential
  end % j
end % k2
stokesDLPtar = stokesDLPtar/pi;
% double-layer potential due to vesicles indexed over K1 
% evaluated at arbitrary points

if nargin == 3
  stokesDLP = zeros(2*vesicle.N,nv);
  for k1 = 1:nv
    K = [(1:k1-1) (k1+1:nv)];
    for j=1:vesicle.N
      diffxy = [X(j,k1) - X(1:vesicle.N,K) ; ...
          X(j+vesicle.N,k1) - X(vesicle.N+1:2*vesicle.N,K)];
      dis2 = diffxy(1:vesicle.N,:).^2 + ...
          diffxy(vesicle.N+1:2*vesicle.N,:).^2;
      % difference of source and target location and distance squared

      coeff = (diffxy(1:vesicle.N,:).*normal(1:vesicle.N,K) + ...
        diffxy(vesicle.N+1:2*vesicle.N,:).*...
        normal(vesicle.N+1:2*vesicle.N,K))./dis2.^2 .* ...
        (diffxy(1:vesicle.N,:).*den(1:vesicle.N,K) + ...
        diffxy(vesicle.N+1:2*vesicle.N,:).*...
        den(vesicle.N+1:2*vesicle.N,K));
      % \frac{(r \dot n)(r \dot density)}{\rho^{4}} term

      stokesDLP(j,k1) = stokesDLP(j,k1) + ...
          sum(sum(coeff.*diffxy(1:vesicle.N,:)));
      stokesDLP(j+vesicle.N,k1) = stokesDLP(j+vesicle.N,k1) + ...
          sum(sum(coeff.*diffxy(vesicle.N+1:2*vesicle.N,:)));
      % double-layer potential for Stokes
    end
  end

  stokesDLP = stokesDLP/pi;
  % 1/pi is the coefficient in front of the double-layer potential
else
  stokesDLP = [];
end
% double-layer potential due to all vesicles except oneself

end % exactStokesDL


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [stokesSLP,stokesSLPtar] = ...
    exactStokesSLfmm(o,vesicle,f,Xtar,K)
% [stokesSLP,stokeSLPtar] = exactStokesSLfmm(vesicle,f,Xtar,K) uses 
% the FMM to compute the single-layer potential due to all vesicles
% except itself vesicle is a class of object capsules and f is the 
% density function NOT scaled by arclength term.  Xtar is a set of 
% points where the single-layer potential due to all vesicles in index 
% set K needs to be evaulated
global fmms

fmms = fmms + 1;
% count the total number of calls to fmm

N = vesicle.N; % number of points per vesicle
nv = vesicle.nv; % number of vesicles
X = vesicle.X; % vesicle positions
oc = curve;
[x,y] = oc.getXY(X); % seperate x and y coordinates

den = f.*[vesicle.sa;vesicle.sa]*2*pi/N;

newFMM = true;
% still want to be able to use the old FMM to check performance

if (nargin == 5)
  stokesSLP = [];
else
  [f1,f2] = oc.getXY(den);
  % need to multiply by arclength term.  Seperate it into
  % x and y coordinate


  if ~newFMM
    potTot = o.fmmStokes(x(:),y(:),f1(:),f2(:));
    % single-layer potential due to all vesicles
    stokesSLP = zeros(2*N,nv); % initialize
    for k = 1:nv
      is = (k-1)*N+1;
      ie = k*N;
      stokesSLP(1:N,k) = potTot(is:ie);
      is = is + N*nv;
      ie = ie + N*nv;
      stokesSLP(N+1:2*N,k) = potTot(is:ie);
    end
    % Wrap the output of the FMM into the usual 
    % [[x1;y1] [x2;y2] ...] format
  end

  if newFMM
    [u,v] = stokesSLPfmm(f1(:),f2(:),x(:),y(:));
    stokesSLP = zeros(2*N,nv); % initialize
    for k = 1:nv
      is = (k-1)*N+1;
      ie = k*N;
      stokesSLP(1:N,k) = u(is:ie);
      stokesSLP(N+1:2*N,k) = v(is:ie);
    end
    % Wrap the output of the FMM into the usual 
    % [[x1;y1] [x2;y2] ...] format
  end


  for k = 1:nv
    if ~newFMM
      pot = o.fmmStokes(x(:,k),y(:,k),f1(:,k),f2(:,k));
      stokesSLP(:,k) = stokesSLP(:,k) - pot;
    end
    if newFMM
      [u,v] = stokesSLPfmm(f1(:,k),f2(:,k),x(:,k),y(:,k));
      stokesSLP(:,k) = stokesSLP(:,k) - [u;v];
    end
  end
  % Subtract potential due to each vesicle on its own.  Nothing
  % to change here for potential at Xtar
end

if nargin == 3
  stokesSLPtar = [];
else
  [x,y] = oc.getXY(X(:,K)); 
  % seperate x and y coordinates at vesicles indexed by K
  [Ntar,ncol] = size(Xtar);
  Ntar = Ntar/2;
  x2 = Xtar(1:Ntar,:);
  x = [x(:);x2(:)];
  y2 = Xtar(Ntar+1:2*Ntar,:);
  y = [y(:);y2(:)];
  % Stack the x and y coordinates of the target points
  [f1,f2] = oc.getXY(den(:,K));
  % seperate x and y coordinates at vesicles indexed by K
  f1 = [f1(:);zeros(Ntar*ncol,1)];
  f2 = [f2(:);zeros(Ntar*ncol,1)];
  % pad density function with zeros so that Xtar doesn't
  % affect the single-layer potential

  if ~newFMM
    potTot = o.fmmStokes(x,y,f1,f2);
    % evalaute single-layer potential using fmm
    stokesSLPtar = zeros(2*Ntar,ncol); % initialize
    for k = 1:ncol
      is = N*numel(K) + (k-1)*Ntar+1;
      ie = is + Ntar - 1;
      stokesSLPtar(1:Ntar,k) = potTot(is:ie);
      is = is + numel(K)*N + Ntar*ncol;
      ie = ie + numel(K)*N + Ntar*ncol;
      stokesSLPtar(Ntar+1:2*Ntar,k) = potTot(is:ie);
    end
    % Wrap the output of the FMM in the usual format
    % for the target points
  end


  if newFMM
    [u,v] = stokesSLPfmm(f1,f2,x,y);
    stokesSLPtar = zeros(2*Ntar,ncol); % initialize
    for k = 1:ncol
      is = N*numel(K) + (k-1)*Ntar+1;
      ie = is + Ntar - 1;
      stokesSLPtar(1:Ntar,k) = u(is:ie);
      stokesSLPtar(Ntar+1:2*Ntar,k) = v(is:ie);
    end
    % Wrap the output of the FMM in the usual format
    % for the target points
  end
end

end % exactStokesSLfmm




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [stokesDLP,stokesDLPtar] = ...
    exactStokesDLfmm(o,vesicle,f,Xtar,K)
% [stokesDLP,stokeDLPtar] = exactStokesDLfmm(vesicle,f,Xtar,K) uses 
% the FMM to compute the double-layer potential due to all vesicles
% except itself vesicle is a class of object capsules and f is the 
% density function NOT scaled by arclength term.  Xtar is a set of 
% points where the single-layer potential due to all vesicles in index 
% set K needs to be evaulated
global fmms

fmms = fmms + 1;
% count the total number of calls to fmm

N = vesicle.N; % number of points per vesicle
nv = vesicle.nv; % number of vesicles
X = vesicle.X; % vesicle positions
oc = curve;
[x,y] = oc.getXY(X); % seperate x and y coordinates
[nx,ny] = oc.getXY(vesicle.normal);
% seperate the x and y coordinates of the normal vector

den = f.*[vesicle.sa;vesicle.sa]*2*pi/N;

if (nargin == 5)
  stokesDLP = [];
else
  [f1,f2] = oc.getXY(den);
  % need to multiply by arclength term.  Seperate it into
  % x and y coordinate

  [u,v] = stokesDLPfmm(f1(:),f2(:),x(:),y(:),nx(:),ny(:));
  stokesDLP = zeros(2*N,nv); % initialize
  for k = 1:nv
    is = (k-1)*N+1;
    ie = k*N;
    stokesDLP(1:N,k) = u(is:ie);
    stokesDLP(N+1:2*N,k) = v(is:ie);
  end
  % Wrap the output of the FMM into the usual 
  % [[x1;y1] [x2;y2] ...] format


  for k = 1:nv
    [u,v] = stokesDLPfmm(f1(:,k),f2(:,k),x(:,k),y(:,k),...
        nx(:,k),ny(:,k));
    stokesDLP(:,k) = stokesDLP(:,k) - [u;v];
  end
  % Subtract potential due to each vesicle on its own.  Nothing
  % to change here for potential at Xtar
end

if nargin == 3
  stokesDLPtar = [];
else
  [x,y] = oc.getXY(X(:,K)); 
  % seperate x and y coordinates at vesicles indexed by K
  [nx,ny] = oc.getXY(vesicle.normal(:,K));
  [Ntar,ncol] = size(Xtar);
  Ntar = Ntar/2;
  x2 = Xtar(1:Ntar,:);
  x = [x(:);x2(:)];
  y2 = Xtar(Ntar+1:2*Ntar,:);
  y = [y(:);y2(:)];
  % Stack the x and y coordinates of the target points
  [f1,f2] = oc.getXY(den(:,K));
  % seperate x and y coordinates at vesicles indexed by K
  f1 = [f1(:);zeros(Ntar*ncol,1)];
  f2 = [f2(:);zeros(Ntar*ncol,1)];
  % pad density function with zeros so that Xtar doesn't
  % affect the double-layer potential
  nx = [nx(:);zeros(Ntar*ncol,1)];
  ny = [ny(:);zeros(Ntar*ncol,1)];
  % pad the normal vector with zeros so that Xtar doesn't
  % affect the double-layer potential

  [u,v] = stokesDLPfmm(f1,f2,x,y,nx,ny);
  stokesDLPtar = zeros(2*Ntar,ncol); % initialize
  for k = 1:ncol
    is = N*numel(K) + (k-1)*Ntar+1;
    ie = is + Ntar - 1;
    stokesDLPtar(1:Ntar,k) = u(is:ie);
    stokesDLPtar(Ntar+1:2*Ntar,k) = v(is:ie);
  end
  % Wrap the output of the FMM in the usual format
  % for the target points
end


end % exactStokesDLfmm




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [laplaceDLP,laplaceDLPtar] = ...
    exactLaplaceDLfmm(o,vesicle,f,Xtar,K)
% [laplaceDLP,laplaceDLPtar] = exactLaplaceDLfmm(vesicle,f) uses 
% the FMM to compute the double-layer laplace potential due to all 
% vesicles except itself vesicle is a class of object capsules and 
% f is the density function NOT scaled by arclength term.  Xtar is 
% a set of points where the single-layer potential due to all 
% vesicles in index set K

N = vesicle.N; % number of points per vesicle
nv = vesicle.nv; % number of vesicles
X = vesicle.X; % vesicle positions
normal = vesicle.normal; % vesicle's outward normal
oc = curve;
[x,y] = oc.getXY(X); % seperate x and y coordinates
[nx,ny] = oc.getXY(normal); % seperate x and y coordinates of normal

[f,~] = oc.getXY(f);
den = f.*vesicle.sa*2*pi/N;

if (nargin == 5)
  laplaceDLP = [];
else
  [~,rfield,~] = fmm_laplace(nx(:).*den(:),x(:),y(:),1);
  [~,~,cfield] = fmm_laplace(ny(:).*den(:),x(:),y(:),1);
  % double-layer potential due to all vesicles
  potTot = rfield + cfield;

  laplaceDLP = zeros(N,nv); % initialize
  for k = 1:nv
    is = (k-1)*N+1;
    ie = k*N;
    laplaceDLP(1:N,k) = potTot(is:ie);
  end
  % Wrap the output of the FMM into the usual 
  % [[x1;y1] [x2;y2] ...] format

  for k = 1:nv
    [~,rfield,~] = fmm_laplace(nx(:,k).*den(:,k),x(:,k),y(:,k),1);
    [~,~,cfield] = fmm_laplace(ny(:,k).*den(:,k),x(:,k),y(:,k),1);
    laplaceDLP(:,k) = laplaceDLP(:,k) - rfield - cfield;
  end
  % Subtract potential due to each vesicle on its own.  Nothing
  % to change here for potential at Xtar
end
laplaceDLP = -[laplaceDLP;zeros(vesicle.N,vesicle.nv)]/(2*pi);
% pad with zeros so that it is compatible with Stokes near-singular
% integration

if nargin == 3
  laplaceDLPtar = [];
else
  x = x(:,K);
  y = y(:,K);
  nx = nx(:,K);
  ny = ny(:,K);
  den = den(:,K);
  % only care about points on vesicles indexed by K
  [Ntar,ncol] = size(Xtar);
  Ntar = Ntar/2;
  x2 = Xtar(1:Ntar,:);
  x = [x(:);x2(:)];
  y2 = Xtar(Ntar+1:2*Ntar,:);
  y = [y(:);y2(:)];
  % Stack the x and y coordinates of the target points
  den = [den(:);zeros(Ntar*ncol,1)];
  nx = [nx(:);zeros(Ntar*ncol,1)];
  ny = [ny(:);zeros(Ntar*ncol,1)];
  % pad density function and normal with zeros so that Xtar doesn't
  % affect the single-layer potential
  [~,rfield,~] = fmm_laplace(nx.*den,x(:),y(:),1);
  [~,~,cfield] = fmm_laplace(ny.*den,x(:),y(:),1);
  % evalaute single-layer potential using fmm

  laplaceDLPtar = zeros(2*Ntar,ncol); % initialize
  for k = 1:ncol
    is = N*numel(K) + (k-1)*Ntar+1;
    ie = is + Ntar - 1;
    laplaceDLPtar(1:Ntar,k) = rfield(is:ie)+cfield(is:ie);
  end
  % Wrap the output of the FMM in the usual format
  % for the target points
  laplaceDLPtar = -laplaceDLPtar/(2*pi);
end


end % exactLaplaceDLfmm


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function pot = fmmStokes(o,x,y,f1,f2)
% pot = fmm_stokes(x,y,f1,f2) computes the single-layer potential
% at points (x,y) due to the density function (f1,f2).  Only diagonal
% term is ommited.

[pot1,rfield,cfield] = fmm_laplace(f1,x,y,2);
pot = [-pot1 + x.*rfield;x.*cfield];
% Logarithm term and part of the r \otimes r term

[pot1,rfield,cfield] = fmm_laplace(f2,x,y,2);
pot = pot + [y.*rfield;-pot1 + y.*cfield];
% Logarithm term and part of the r \otimes r term

[~,rfield,cfield] = fmm_laplace(x.*f1 + y.*f2,x,y,1);
pot = pot - [rfield;cfield];
% Other part of r \otimes r term

pot = 0.25*pot/pi;
% single-layer potential has 1/4/pi term in front


end % fmmStokes




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [laplaceDLP,laplaceDLPtar] = ...
    exactLaplaceDL(o,vesicle,f,Xtar,K1)
% pot = exactLaplaceDL(vesicle,f,Xtar,K1) 
% computes the double-layer laplace potential due to f around all
% vesicles except itself.  Also can pass a set of target points Xtar
% and a collection of vesicles K1 and the double-layer potential due 
% to vesicles in K1 will be evaluated at Xtar.
% Everything but Xtar is in the 2*N x nv format
% Xtar is in the 2*Ntar x ncol format

nv = vesicle.nv; % number of vesicles
X = vesicle.X; % Vesicle positions
normal = vesicle.normal; % Outward normal
sa = vesicle.sa; % Jacobian

if nargin == 5
  Ntar = size(Xtar,1)/2;
  ncol = size(Xtar,2);
  laplaceDLPtar = zeros(2*Ntar,ncol);
else
  K1 = [];
  laplaceDLPtar = [];
  ncol = 0;
  Ntar = 0;
  % if nargin ~= 5, user does not need the layer potential at arbitrary
  % points
end

den = f.*[sa;sa]*2*pi/vesicle.N;
% multiply by arclength term

for k2 = 1:ncol % loop over columns of target points
  for j = 1:Ntar % loop over rows of target points
    dis2 = (X(1:vesicle.N,K1) - Xtar(j,k2)).^2 + ... 
        (X(vesicle.N+1:2*vesicle.N,K1) - Xtar(j+Ntar,k2)).^2;
    diffxy = [X(1:vesicle.N,K1) - Xtar(j,k2) ; ...
        X(vesicle.N+1:2*vesicle.N,K1) - Xtar(j+Ntar,k2)];
    % distance squared and difference of source and target location

    coeff = (diffxy(1:vesicle.N,:).*normal(1:vesicle.N,K1) + ...
      diffxy(vesicle.N+1:2*vesicle.N,:).*...
      normal(vesicle.N+1:2*vesicle.N,K1))./dis2;
    % this is the kernel of the double-layer potential for
    % laplace's equation
    val = coeff.*den(1:vesicle.N,K1);
    laplaceDLPtar(j,k2) = sum(val(:));
    val = coeff.*den(vesicle.N+1:2*vesicle.N,K1);
    laplaceDLPtar(j+Ntar,k2) = sum(val(:));
  end % j

end % k2
% Evaluate double-layer potential at arbitrary target points
laplaceDLPtar = 1/(2*pi)*laplaceDLPtar;
% 1/2/pi is the coefficient in front of the double-layer potential

laplaceDLP = zeros(vesicle.N,vesicle.nv); % Initialize to zero
% if we only have one vesicle, vesicles of course can not collide
% Don't need to run this loop in this case
if (nargin == 3 && vesicle.nv > 1)
  for k1 = 1:vesicle.nv % vesicle of targets
    K = [(1:k1-1) (k1+1:vesicle.nv)];
    % Loop over all vesicles except k1

    for j=1:vesicle.N
      dis2 = (X(1:vesicle.N,K) - X(j,k1)).^2 + ...
          (X(vesicle.N+1:2*vesicle.N,K) - X(j+vesicle.N,k1)).^2;
      diffxy = [X(1:vesicle.N,K) - X(j,k1); ...
          X(vesicle.N+1:2*vesicle.N,K) - X(j+vesicle.N,k1)];
      % distance squared and difference of source and target location

      coeff = (diffxy(1:vesicle.N,:).*normal(1:vesicle.N,K) + ...
        diffxy(vesicle.N+1:2*vesicle.N,:).*...
        normal(vesicle.N+1:2*vesicle.N,K))./dis2;
      % this is the kernel of the double-layer potential for
      % laplace's equation
      val = coeff.*den(1:vesicle.N,K);
      laplaceDLP(j,k1) = sum(val(:));
    end % j
  end % k1
  % Evaluate double-layer potential at vesicles but oneself
  laplaceDLP = 1/(2*pi)*laplaceDLP;
  % 1/2/pi is the coefficient in front of the double-layer potential
end % nargin == 3

laplaceDLP = [laplaceDLP;laplaceDLP];

end % exactLaplaceDL




end % methods

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
methods (Static)

function LP = lagrangeInterp
% interpMap = lagrangeInterp builds the Lagrange interpolation
% matrix that takes seven function values equally distributed
% in [0,1] and returns the seven polynomial coefficients

interpMat = zeros(7);
LP(1,1) = 6.48e1;
LP(1,2) = -3.888e2;
LP(1,3) = 9.72e2;
LP(1,4) = -1.296e3;
LP(1,5) = 9.72e2;
LP(1,6) = -3.888e2;
LP(1,7) = 6.48e1;

LP(2,1) = -2.268e2;
LP(2,2) = 1.296e3;
LP(2,3) = -3.078e3;
LP(2,4) = 3.888e3;
LP(2,5) = -2.754e3;
LP(2,6) = 1.0368e3;
LP(2,7) = -1.62e2;

LP(3,1) = 3.15e2;
LP(3,2) = -1.674e3;
LP(3,3) = 3.699e3;
LP(3,4) = -4.356e3;
LP(3,5) = 2.889e3;
LP(3,6) = -1.026e3;
LP(3,7) = 1.53e2;

LP(4,1) = -2.205e2;
LP(4,2) = 1.044e3;
LP(4,3) = -2.0745e3;
LP(4,4) = 2.232e3;
LP(4,5) = -1.3815e3;
LP(4,6) = 4.68e2;
LP(4,7) = -6.75e1;

LP(5,1) = 8.12e1;
LP(5,2) = -3.132e2;
LP(5,3) = 5.265e2;
LP(5,4) = -5.08e2;
LP(5,5) = 2.97e2;
LP(5,6) = -9.72e1;
LP(5,7) = 1.37e1;

LP(6,1) = -1.47e1;
LP(6,2) = 3.6e1;
LP(6,3) = -4.5e1;
LP(6,4) = 4.0e1;
LP(6,5) = -2.25e1;
LP(6,6) = 7.2e0;
LP(6,7) = -1e0;

LP(7,1) = 1e0;
% rest of the coefficients are zero


end % lagrangeInterp


end % methods (Static)





end % classdef
