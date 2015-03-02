classdef curve 
% This class implements that basic calculus on the curve.
% The basic data structure is a matrix X in which the columns 
% represent periodic C^{\infty} closed curves with N points, 
% X(1:n,j) is the x-coordinate of the j_th curve and X(n+1:N,j) 
% is the y-coordinate of the j_th curve; here n=N/2
% X coordinates do not have to be periodic, but the curvature,
% normals, etc that they compute will be garbage.

methods

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [x,y]=getXY(o,X)
% [x,y] = getXY(X) get the [x,y] component of curves X
N = size(X,1)/2;
x = X(1:N,:);
y = X(N+1:end,:);
end % getXY


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function V = setXY(o,x,y)
% V = setXY(x,y) set the [x,y] component of vector V on the curve
N = size(x,1);
V=zeros(2*N,size(x,2));
V(1:N,:) = x;
V(N+1:end,:) = y;
end % setXY


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [Dx,Dy]=getDXY(o,X)
% [Dx,Dy]=getDXY(X), compute the derivatives of each component of X 
% these are the derivatives with respect the parameterization 
% not arclength
[x,y] = o.getXY(X);
N = size(x,1);
nv = size(x,2);
IK = o.modes(N,nv);
Dx = o.diffFT(x,IK);
Dy = o.diffFT(y,IK);
end % getDXY


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [jacobian,tangent,curvature] = diffProp(o,X)
% [jacobian,tangent,curvature] = diffProp(X) returns differential
% properties of the curve each column of the matrix X. Each column of 
% X should be a closed curve defined in plane. The tangent is the 
% normalized tangent.
%
n = size(X,1);
nv = size(X,2);
N = n/2; 
IK = o.modes(N,nv);

% get the x y components
[Dx,Dy] = o.getDXY(X);

jacobian = sqrt( Dx.^2 + Dy.^2 ); 

if nargout>1  % if user requires tangent
  tangent = o.setXY( Dx./jacobian, Dy./jacobian);
end

if nargout>2  % if user requires curvature
  DDx = curve.arcDeriv(Dx,1,ones(N,1),IK);
  DDy = curve.arcDeriv(Dy,1,ones(N,1),IK);
  curvature = (Dx.*DDy - Dy.*DDx)./(jacobian.^3);
end
% curvature of the curve

end % diffProp

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [reducedArea,area,length] = geomProp(o,X)
% [reducedArea area length] = geomProp(X) calculate the length, area 
% and the reduced volume of domains inclose by columns of X. 
% Reduced volume is defined as 4*pi*A/L. 

[x,y] = o.getXY(X);
N = size(x,1);
[Dx,Dy] = o.getDXY(X);

length = sum( sqrt(Dx.^2 + Dy.^2) )*2*pi/N;
area = sum( x.*Dy - y.*Dx)*pi/N;

reducedArea = 4*pi*area./length.^2;

end % geomProp


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function A = sinterpS(o,N,y)
% A = sinterpS(N,y) constructs the interpolation matrix A that maps
% a function defined periodically at N equispaced points to the
% function value at the points y
% The points y are assumed to be the in the 0-2*pi range. 

A = zeros(numel(y),N);
modes = [(0:N/2-1) 0 (-N/2+1:-1)];
for j=1:N
  f = zeros(N,1); f(j) = 1;
  fhat = fft(f)/N;
  for k=1:N
    A(:,j) = A(:,j) + fhat(k)*exp(1i*modes(k)*y);
  end
end
A = real(A);
% Input is real, so interpolant should be real

end % sinterpS


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [X,nv] = initConfig(o,N,varargin)       
% [X,nv] = initConfig(n,varargin) returns N coordinates of boundary
% points.
% X = BOUNDARY(N,OPTIONS) can be used to call for different
% configuration.  Available OPTIONS are
%
%   'nv'          - (followed by) number of vesicles (in 
%                   linear configuration),
%   'circles'     - circle of given radius and center
%                   position, 
%   'square'      - returns a rounded square domain 

options = varargin;

if(any(strcmp(options,'nv')))
  nv = options{find(strcmp(options,'nv'))+1};
else
  nv = 1;
end
% number of components

if(any(strcmp(options,'radii')))
  radii = options{find(strcmp(options,'radii'))+1};
else
  radii = ones(nv,1);
end
% radii of the boudnary curves

if(any(strcmp(options,'center')))
  center = options{find(strcmp(options,'center'))+1};
else
  center = [(0:nv-1);zeros(1,nv)];
end
% centers of the boundary curve

if(any(strcmp(options,'radiiBeans')))
  radiiBeans = options{find(strcmp(options,'radiiBeans'))+1};
else
  radiiBeans = [];
end
% radii of posts cut out to form beans

if(any(strcmp(options,'centerBeans')))
  centerBeans = options{find(strcmp(options,'centerBeans'))+1};
else
  centerBeans = [];
end
% centers of posts cut out to form beans

t = (0:N-1)'*2*pi/N;
% Discretization in parameter space
X = zeros(2*N,nv);

if any(strcmp(options,'circles'))
  for k = 1:nv
    X(:,k) = [center(k,1) + radii(k)*cos(t); ...
              center(k,2) + radii(k)*sin(t)];
  end
  % paramterization of the circular exclusions

elseif any(strcmp(options,'beans'))
  indCutGrains = [360 139 403 195 116 253 323 398 263 364 362 223 306 ...
      349 431 179 1 396 109 315 418 277 295 374 310 183 260 ...
      289 447 185 313 338 143 404 307 128 280 347 324 406 420 316 ...
      287 416 434 370 436 312 432 402 415 384 369 405 408 333 417 ...
      410 373 351 443 371 353 279 321 298 411 422 392 388 231 424 ...
      345 445 452 236 429 217 400 401 439 379 234 380 309 412 461 ...
      361 302 433 464 221 346 435 459 414 419 391 423 383 442 219 ...
      428 437 375 465 305 451 457 441 462 395 440 367 463 430 444 ...
      448 446 427 449 460 409 454 450]; % 421 421];
  indFullGrains = setdiff((1:465),[indCutGrains 421]);
  for k = 1:numel(indFullGrains)
    kk = indFullGrains(k);
    X(:,k) = [center(kk,1) + radii(kk)*cos(t); ...
              center(kk,2) + radii(kk)*sin(t)];
  end

  for k = 1:numel(indCutGrains);
    kk = indCutGrains(k);
    z = intersectingCircles(N,radii(kk),center(kk,1),center(kk,2),...
        radiiBeans(k),centerBeans(k,1),centerBeans(k,2));
    modes = (-N/2:N/2-1)';
    zh = fftshift(fft(z));
    zh(abs(modes) > 20) = 0;
    z = ifft(ifftshift(zh));
    X(:,k+numel(indFullGrains)) = [real(z);imag(z)];
  end

  kk = 421;
  k1 = 126;
  k2 = 127;
  z = appleGeom(N,radii(kk),center(kk,1),center(kk,2),...
      radiiBeans(k1),centerBeans(k1,1),centerBeans(k1,2),...
      radiiBeans(k2),centerBeans(k2,1),centerBeans(k2,2));
  modes = (-N/2:N/2-1)';
  zh = fftshift(fft(z));
  zh(abs(modes) > 20) = 0;
  z = ifft(ifftshift(zh));
  X(:,465) = [real(z);imag(z)];

%  for k = 1:numel(radiiBeans)
%    X(:,k+465) = [centerBeans(k,1) + radiiBeans(k)*cos(t); ...
%              centerBeans(k,2) + radiiBeans(k)*sin(t)];
%  end




  % index of grains/posts that are cut by a grain in
  % radiiBeans/centersBeans.  Order is so that bean(k) cuts through
  % grain indCutGrains(k)

elseif any(strcmp(options,'square'))
  a = 22.51; b = 2.6e0; order = 10;
  % parameters for the boundary
  r = (cos(t).^order + sin(t).^order).^(-1/order);
  x = a*r.*(cos(t))+a; y = b*r.*(sin(t))+b;
  x = x - 7.7;

  X = [x;y];
  % rounded off square.  Increase order ot make it more square

elseif any(strcmp(options,'square12'))
  a = 25.8; b = 2.6e0; order = 10;
  % parameters for the boundary
  r = (cos(t).^order + sin(t).^order).^(-1/order);
  x = a*r.*(cos(t))+a; y = b*r.*(sin(t))+b;
  x = x - 7.8811;

  X = [x;y];
  % rounded off square.  Increase order ot make it more square

elseif any(strcmp(options,'square17'))
  a = 18.9; b = 2.6e0; order = 10;
  % parameters for the boundary
  r = (cos(t).^order + sin(t).^order).^(-1/order);
  x = a*r.*(cos(t))+a; y = b*r.*(sin(t))+b;
  x = x - 5.8;

  X = [x;y];
  % rounded off square.  Increase order ot make it more square

elseif any(strcmp(options,'square20'))
  a = 25.3; b = 2.6e0; order = 10;
  % parameters for the boundary
  r = (cos(t).^order + sin(t).^order).^(-1/order);
  x = a*r.*(cos(t))+a; y = b*r.*(sin(t))+b;
  x = x - 7.75336;

  X = [x;y];
  % rounded off square.  Increase order ot make it more square
elseif any(strcmp(options,'square40'))
  a = 18.9; b = 2.6e0; order = 10;
  % parameters for the boundary
  r = (cos(t).^order + sin(t).^order).^(-1/order);
  x = a*r.*(cos(t))+a; y = b*r.*(sin(t))+b;
  x = x - 5.8;

  X = [x;y];
  % rounded off square.  Increase order ot make it more square

elseif any(strcmp(options,'square45'))
  a = 18.9; b = 2.6e0; order = 10;
  % parameters for the boundary
  r = (cos(t).^order + sin(t).^order).^(-1/order);
  x = a*r.*(cos(t))+a; y = b*r.*(sin(t))+b;
  x = x - 5.8;

  X = [x;y];
  % rounded off square.  Increase order ot make it more square

elseif any(strcmp(options,'square46'))
  a = 14.2; b = 2.6e0; order = 10;
  % parameters for the boundary
  r = (cos(t).^order + sin(t).^order).^(-1/order);
  x = a*r.*(cos(t))+a; y = b*r.*(sin(t))+b;
  x = x - 4.5;

  X = [x;y];
  % rounded off square.  Increase order ot make it more square

end


end % initConfig


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function indOut = collision(o,Xtar,geom,Near,fmm)
% indOut = collision(o,Xtar,geom,Near,fmm)

[x,y] = o.getXY(geom.X);
[nx,ny] = o.getXY(geom.normal);

f = [ones(geom.N,geom.nv);zeros(geom.N,geom.nv)];
% Density function is constant.  Pad second half of it with zero
op = poten(geom,fmm);
% load object for doing near-singular integration and evaluating
% laplace double-layer potential

if ~fmm
  kernel = @op.exactLaplaceDL;
else
  kernel = @op.exactLaplaceDL;
%  kernel = @op.exactLaplaceDLfmm;
% TODO: Need to put this in
end
% kernel for laplace's double layer potential.  Only difference
% is if FMM is used or not

Fdlp = op.nearSingInt(geom,f,@(X,f) 1.5*f,...
       Near,kernel,Xtar,false,'outer');
% take inner since we are already accounting for the jump by assiging 0
% since the density function is 1.
Fdlp = Fdlp(1:end/2);

thresh = 9.999e-1;
indOut = find(abs(Fdlp) < thresh);

end % collision






end % methods


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
methods (Static)

function df = arcDeriv(f,m,sa,IK,trash)
% df = arcDeriv(f,m,s,IK) is the arclength derivative of order m.
% f is a matrix of scalar functions (each function is a column)
% f is assumed to have an arbitrary parametrization
% sa = d s/ d a, where a is the aribtrary parameterization
% IK is the fourier modes which is saved and used to accelerate 
% this routine

if nargin<3, m=1; end;
if nargin<4, sa=ones(size(f,1),1); end;

col = size(f,2);
if col == 2
  sa = [sa sa];
elseif col == 3
  sa = [sa sa sa];
elseif col > 3
  sa = repmat(sa,1,col);
end
% This is much faster than always using repmat


%if nargin == 5
%  frac = 1;
%  N = size(IK,1);
%  op = poten(N);
%  IK(abs(IK) > frac*N/2) = 0;
%  norm(o.TimeMatVec(Xn,vesicle,walls,op)-rhs,inf)/norm(rhs,inf)
%  pause
%end

df = f;
for j=1:m
  df = sa.*real(ifft(IK.*fft(df)));
end
% compute the mth order derivative

end % arcDeriv


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function Deriv = D1(N)
% Deriv = D1(N) constructs a N by N fourier differentiation matrix
[FF,FFI] = fft1.fourierInt(N);
Deriv = FFI * diag(1i*([0 -N/2+1:N/2-1])) * FF;
Deriv = real(Deriv);

end % D1


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function df = diffFT(f,IK)
% df = diffFT(f,IK) Computes the first derivative of an array of 
% periodic functions f using fourier transform. The f(:,1) is the 
% first function, f(:,2) is the second function, etc.
% IK is used to speed up the code.  It is the index of the fourier
% modes so that fft and ifft can be used
% 

df = real(ifft(IK.*fft(f)));

end % end diffFT


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function IK = modes(N,nv)
% IK = modes(N) builds the order of the fourier modes required for using
% fft and ifft to do spectral differentiation

IK = 1i*[(0:N/2-1) 0 (-N/2+1:-1)]';
% diagonal term for Fourier differentiation with the -N/2 mode
% zeroed to avoid Nyquist frequency

if nv == 2
  IK = [IK IK];
elseif nv == 3
  IK = [IK IK IK];
elseif nv > 3
  IK = repmat(IK,1,nv);
end

end % modes


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [FF,FFI] = fourierInt(N)
% [FF,FFI] = fourierInt(N) returns a matrix that take in the point
% values of a function in [0,2*pi) and returns the fourier coefficients % (FF) and a matrix that takes in the fourier coefficients and returns
% the function values (FFI)

theta = (0:N-1)'*2*pi/N;
%modes = [0;(-N/2+1:N/2-1)'];
modes = [-N/2;(-N/2+1:N/2-1)'];

FF = zeros(N);
for i=1:N
  FF(:,i) = exp(-1i*theta(i)*modes);
end
FF = FF/N;
% FF takes function values and returns the Fourier coefficients

if (nargout > 1)
  FFI = zeros(N);
  for i=1:N
    FFI(:,i) = exp(1i*modes(i)*theta);
  end
  % FFI takes the Fourier coefficients and returns function values.
else
  FFI = [];
end

end % fourierInt

end % methods (Static)


end % classdef


