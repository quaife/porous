function zi = interp2FAST(x,y,z,xi,yi)
%INTERP2 2-D interpolation (table lookup).
%   ZI = INTERP2(X,Y,Z,XI,YI) interpolates to find ZI, the values of the
%   underlying 2-D function Z at the points in matrices XI and YI.
%   Matrices X and Y specify the points at which the data Z is given.
%
%   XI can be a row vector, in which case it specifies a matrix with
%   constant columns. Similarly, YI can be a column vector and it 
%   specifies a matrix with constant rows. 
%
%   ZI = INTERP2(Z,XI,YI) assumes X=1:N and Y=1:M where [M,N]=SIZE(Z).
%   ZI = INTERP2(Z,NTIMES) expands Z by interleaving interpolates between
%   every element, working recursively for NTIMES.  INTERP2(Z) is the
%   same as INTERP2(Z,1).
%
%   ZI = INTERP2(...,METHOD) specifies alternate methods.  The default
%   is linear interpolation.  Available methods are:
%
%     'nearest' - nearest neighbor interpolation
%     'linear'  - bilinear interpolation
%     'spline'  - spline interpolation
%     'cubic'   - bicubic interpolation as long as the data is
%                 uniformly spaced, otherwise the same as 'spline'
%
%   For faster interpolation when X and Y are equally spaced and monotonic,
%   use the syntax ZI = INTERP2(...,*METHOD).
%
%   ZI = INTERP2(...,METHOD,EXTRAPVAL) specificies a method and a scalar 
%   value for ZI outside of the domain created by X and Y.  Thus, ZI will
%   equal EXTRAPVAL for any value of YI or XI which is not spanned by Y 
%   or X respectively. A method must be specified for EXTRAPVAL to be used,
%   the default method is 'linear'.
%
%   All the interpolation methods require that X and Y be monotonic and
%   plaid (as if they were created using MESHGRID).  If you provide two
%   monotonic vectors, interp2 changes them to a plaid internally. 
%   X and Y can be non-uniformly spaced.
%
%   For example, to generate a coarse approximation of PEAKS and
%   interpolate over a finer mesh:
%       [x,y,z] = peaks(10); [xi,yi] = meshgrid(-3:.1:3,-3:.1:3);
%       zi = interp2(x,y,z,xi,yi); mesh(xi,yi,zi)
%
%   Class support for inputs X, Y, Z, XI, YI:  
%      float: double, single
%
%   See also INTERP1, INTERP3, INTERPN, MESHGRID, TriScatteredInterp.

%   Copyright 1984-2010 The MathWorks, Inc.
%   $Revision: 5.33.4.23 $  $Date: 2010/11/17 11:29:29 $


%narg = nargin-1;
%method = [varargin{end} '    ']; % Protect against short string.
ExtrapVal = nan; % setting default ExtrapVal as NAN


%[msg,x,y,z,xi,yi] = xyzchk(varargin{1:5});
%x = varargin{1};
%y = varargin{2};
%z = varargin{3};
%xi = varargin{4};
%yi = varargin{5};

zi = cubic(ExtrapVal,x,y,z,xi,yi);

%end


%------------------------------------------------------
function F = cubic(ExtrapVal,arg1,arg2,arg3,arg4,arg5)
%CUBIC 2-D bicubic data interpolation.
%   CUBIC(...) is the same as LINEAR(....) except that it uses
%   bicubic interpolation.
%
%   This function needs about 7-8 times SIZE(XI) memory to be available.
%
%   See also LINEAR.

%   Based on "Cubic Convolution Interpolation for Digital Image
%   Processing", Robert G. Keys, IEEE Trans. on Acoustics, Speech, and
%   Signal Processing, Vol. 29, No. 6, Dec. 1981, pp. 1153-1160.

[nrows,ncols] = size(arg3);
%mx = numel(arg1); my = numel(arg2);
s = 1 + (arg4-arg1(1))/(arg1(end)-arg1(1))*(ncols-1);
t = 1 + (arg5-arg2(1))/(arg2(end)-arg2(1))*(nrows-1);

% Check for out of range values of s and set to 1
%sout = find((s<1)|(s>ncols));

% Check for out of range values of t and set to 1
%tout = find((t<1)|(t>nrows));

% Matrix element indexing
ndx = floor(t)+floor(s-1)*(nrows+2);

% Compute intepolation parameters, check for boundary value.
%d = find(s==ncols);
s(:) = (s - floor(s));

% Compute intepolation parameters, check for boundary value.
%d = find(t==nrows);
t(:) = (t - floor(t));

% Expand z so interpolation is valid at the boundaries.
zz = zeros(size(arg3)+2);
zz(1,2:ncols+1) = 3*arg3(1,:)-3*arg3(2,:)+arg3(3,:);
zz(2:nrows+1,2:ncols+1) = arg3;
zz(nrows+2,2:ncols+1) = 3*arg3(nrows,:)-...
    3*arg3(nrows-1,:)+arg3(nrows-2,:);
zz(:,1) = 3*zz(:,2)-3*zz(:,3)+zz(:,4);
zz(:,ncols+2) = 3*zz(:,ncols+1)-3*zz(:,ncols)+zz(:,ncols-1);
nrows = nrows+2; %also ncols = ncols+2;

% Now interpolate using computationally efficient algorithm.
t0 = ((2-t).*t-1).*t;
t1 = (3*t-5).*t.*t+2;
t2 = ((4-3*t).*t+1).*t;
t(:) = (t-1).*t.*t;
F = ( zz(ndx).*t0 + zz(ndx+1).*t1 + zz(ndx+2).*t2 + ...
    zz(ndx+3).*t ).* (((2-s).*s-1).*s);
ndx(:) = ndx + nrows;
F(:) = F + ( zz(ndx).*t0 + zz(ndx+1).*t1 + zz(ndx+2).*t2 + ...
    zz(ndx+3).*t ).* ((3*s-5).*s.*s+2);
ndx(:) = ndx + nrows;
F(:) = F + ( zz(ndx).*t0 + zz(ndx+1).*t1 + zz(ndx+2).*t2 + ...
    zz(ndx+3).*t ).* (((4-3*s).*s+1).*s);
ndx(:) = ndx + nrows;
F(:) = F + ( zz(ndx).*t0 + zz(ndx+1).*t1 + zz(ndx+2).*t2 + ....
    zz(ndx+3).*t ).* ((s-1).*s.*s);
F(:) = 0.25*F;
% Now set out of range values to ExtrapVal.


%----------------------------------------------------------
function F = spline2(varargin)
%2-D spline interpolation

% Determine abscissa vectors
varargin{1} = varargin{1}(1,:);
varargin{2} = varargin{2}(:,1).';

%
% Check for plaid data.
%
xi = varargin{4}; yi = varargin{5};
xxi = xi(1,:); yyi = yi(:,1);

if ~isequal(repmat(xxi,size(xi,1),1),xi) || ...
        ~isequal(repmat(yyi,1,size(yi,2)),yi)
    F = splncore(varargin(2:-1:1),varargin{3},varargin(5:-1:4));
else
    F = splncore(varargin(2:-1:1),varargin{3},{yyi(:).' xxi},'gridded');
end

ExtrapVal = varargin{6};
% Set out-of-range values to ExtrapVal
if isnumeric(ExtrapVal)
    d = xi < min(varargin{1}) | xi > max(varargin{1}) | ...
        yi < min(varargin{2}) | yi > max(varargin{2});
    F(d) = ExtrapVal;
end
