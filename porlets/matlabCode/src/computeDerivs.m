function [u_x,u_y,v_x,v_y] = computeDerivs(eX,eY,u,v);
% find gradient of velocity field using finite differences.  Want to use
% one-sided derivatives on the boundary and when the stencil is cut by a
% pore

order = 4;
load radii.dat;
load centers.dat;
[ny,nx] = size(eX);

u_x = zeros(ny,nx);
u_y = u_x; v_x = u_x; v_y = u_x;
% allocate space for derivatives
dx = eX(1,2) - eX(1,1);
dy = eY(2,1) - eY(1,1);
% spacing in x and y direction

if order == 1
  % interior points
  u_x(:,1:end-1) = (u(:,2:end) - u(:,1:end-1))/dx;
  u_y(1:end-1,:) = (u(2:end,:) - u(1:end-1,:))/dy;
  v_x(:,1:end-1) = (v(:,2:end) - v(:,1:end-1))/dx;
  v_y(1:end-1,:) = (v(2:end,:) - v(1:end-1,:))/dy;

  % boundary points
  u_x(:,end) = (u(:,end) - u(:,end-1))/dx;
  u_y(end,:) = (u(end,:) - u(end-1,:))/dy;
  v_x(:,end) = (v(:,end) - v(:,end-1))/dx;
  v_y(end,:) = (v(end,:) - v(end-1,:))/dy;

  % Fix points that border a pore (ie. cut cells)
  % Allocate space for finding points whose stencil may cut through a
  % pore
  extPtsX = [];
  extPtsY = [];
  for k = 1:numel(radii)
    dist2 = (eX - centers(k,1)).^2 + (eY - centers(k,2)).^2;
    intPts = find(dist2 > (radii(k)-max(dx,dy))^2 & dist2 <= radii(k)^2);
    u_x(intPts) = 0;
    u_y(intPts) = 0;
    v_x(intPts) = 0;
    v_y(intPts) = 0;
    % Assign zero in the interior

    [newPtsY,newPtsX] = find(dist2 < ...
        (radii(k)+max(dx,dy))^2 & dist2 > radii(k)^2);
    extPtsX = [extPtsX; newPtsX]; 
    extPtsY = [extPtsY; newPtsY]; 
  end

  for k = 1:numel(extPtsX)
    j = extPtsY(k);
    i = extPtsX(k);
    % points whose right neighbor to the right is inside a pore
    if any(((eX(j,i+1) - centers(:,1)).^2 + ...
          (eY(j,i) - centers(:,2)).^2) < radii.^2)
      u_x(j,i) = -(u(j,i-1) - u(j,i))/dx;
      v_x(j,i) = -(v(j,i-1) - v(j,i))/dx;
    end
    % points whose left neighbor to the right is inside a pore
    if any(((eX(j,i-1) - centers(:,1)).^2 + ...
          (eY(j,i) - centers(:,2)).^2) < radii.^2)
      u_x(j,i) = (u(j,i+1) - u(j,i))/dx;
      v_x(j,i) = (v(j,i+1) - v(j,i))/dx;
    end
    % points whose upper neighbor to the right is inside a pore
    if any(((eX(j,i) - centers(:,1)).^2 + ...
          (eY(j+1,i) - centers(:,2)).^2) < radii.^2)
      u_y(j,i) = -(u(j-1,i) - u(j,i))/dx;
      v_y(j,i) = -(v(j-1,i) - v(j,i))/dx;
    end
    % points whose lower neighbor to the right is inside a pore
    if any(((eX(j,i) - centers(:,1)).^2 + ...
          (eY(j-1,i) - centers(:,2)).^2) < radii.^2)
      u_y(j,i) = (u(j+1,i) - u(j,i))/dx;
      v_y(j,i) = (v(j+1,i) - v(j,i))/dx;
    end
  end
end % order == 1

if order == 2
  % Second-order at interior nodes
  u_x(:,2:end-1) = (u(:,3:end) - u(:,1:end-2))/2/dx;
  u_y(2:end-1,:) = (u(3:end,:) - u(1:end-2,:))/2/dy;
  v_x(:,2:end-1) = (v(:,3:end) - v(:,1:end-2))/2/dx;
  v_y(2:end-1,:) = (v(3:end,:) - v(1:end-2,:))/2/dy;
  % second-order finite difference for every point except the final
  % columns/rows

  % Second-order at boundary nodes
  u_x(:,1) = (-0.5*u(:,3) + 2*u(:,2) - 1.5*u(:,1))/dx;
  v_x(:,1) = (-0.5*v(:,3) + 2*v(:,2) - 1.5*v(:,1))/dx;
  u_y(1,:) = (-0.5*u(3,:) + 2*u(2,:) - 1.5*u(1,:))/dy;
  v_y(1,:) = (-0.5*v(3,:) + 2*v(2,:) - 1.5*v(1,:))/dy;
  u_x(:,end) = -(-0.5*u(:,end-2) + 2*u(:,end-1) - 1.5*u(:,end))/dx;
  v_x(:,end) = -(-0.5*v(:,end-2) + 2*v(:,end-1) - 1.5*v(:,end))/dx;
  u_y(end,:) = -(-0.5*u(end-2,:) + 2*u(end-1,:) - 1.5*u(end,:))/dy;
  v_y(end,:) = -(-0.5*v(end-2,:) + 2*v(end-1,:) - 1.5*v(end,:))/dy;


  % Fix points that border a pore (ie. cut cells)
  % Space for finding points whose stencil may cut through a pore
  extPtsX = [];
  extPtsY = [];
  for k = 1:numel(radii)
%    disp(numel(radii) - k + 1)
    dist2 = (eX - centers(k,1)).^2 + (eY - centers(k,2)).^2;
    intPts = find(dist2 > (radii(k)-max(dx,dy))^2 & dist2 <= radii(k)^2);
    u_x(intPts) = 0;
    u_y(intPts) = 0;
    v_x(intPts) = 0;
    v_y(intPts) = 0;
    % Assign zero in the interior

    [newPtsY,newPtsX] = find(dist2 < ...
        (radii(k)+max(dx,dy))^2 & dist2 > radii(k)^2);
    extPtsX = [extPtsX; newPtsX]; 
    extPtsY = [extPtsY; newPtsY]; 
  end

  for k = 1:numel(extPtsX)
    j = extPtsY(k);
    i = extPtsX(k);
    % points whose right neighbor to the right is inside a pore
    if any(((eX(j,i+1) - centers(:,1)).^2 + ...
          (eY(j,i) - centers(:,2)).^2) < radii.^2)
      u_x(j,i) = -(-0.5*u(j,i-2) + 2*u(j,i-1) - 1.5*u(j,i))/dx;
      v_x(j,i) = -(-0.5*v(j,i-2) + 2*v(j,i-1) - 1.5*v(j,i))/dx;
    end
    % points whose left neighbor to the right is inside a pore
    if any(((eX(j,i-1) - centers(:,1)).^2 + ...
          (eY(j,i) - centers(:,2)).^2) < radii.^2)
      u_x(j,i) = (-0.5*u(j,i+2) + 2*u(j,i+1) - 1.5*u(j,i))/dx;
      v_x(j,i) = (-0.5*v(j,i+2) + 2*v(j,i+1) - 1.5*v(j,i))/dx;
    end
    % points whose upper neighbor to the right is inside a pore
    if any(((eX(j,i) - centers(:,1)).^2 + ...
          (eY(j+1,i) - centers(:,2)).^2) < radii.^2)
      u_y(j,i) = -(-0.5*u(j-2,i) + 2*u(j-1,i) - 1.5*u(j,i))/dy;
      v_y(j,i) = -(-0.5*v(j-2,i) + 2*v(j-1,i) - 1.5*v(j,i))/dy;
    end
    % points whose lower neighbor to the right is inside a pore
    if any(((eX(j,i) - centers(:,1)).^2 + ...
          (eY(j-1,i) - centers(:,2)).^2) < radii.^2)
      u_y(j,i) = (-0.5*u(j+2,i) + 2*u(j+1,i) - 1.5*u(j,i))/dy;
      v_y(j,i) = (-0.5*v(j+2,i) + 2*v(j+1,i) - 1.5*v(j,i))/dy;
    end
  end

end % order == 2

if order == 4
  % Fourth-order at interior nodes
  u_x(:,3:end-2) = (-1/12*u(:,5:end) + 2/3*u(:,4:end-1) - ...
      2/3*u(:,2:end-3) + 1/12*u(:,1:end-4))/dx;
  u_y(3:end-2,:) = (-1/12*u(5:end,:) + 2/3*u(4:end-1,:) - ...
      2/3*u(2:end-3,:) + 1/12*u(1:end-4,:))/dy;
  v_x(:,3:end-2) = (-1/12*v(:,5:end) + 2/3*v(:,4:end-1) - ...
      2/3*v(:,2:end-3) + 1/12*v(:,1:end-4))/dx;
  v_y(3:end-2,:) = (-1/12*v(5:end,:) + 2/3*v(4:end-1,:) - ...
      2/3*v(2:end-3,:) + 1/12*v(1:end-4,:))/dy;
  % fourth-order finite difference for every point except the final
  % columns/rows

%  % Fourth-order at boundary nodes
  u_x(:,1) = (-1/4*u(:,5) + 4/3*u(:,4) - ...
      3*u(:,3) + 4*u(:,2) - 25/12*u(:,1))/dx;
  u_x(:,2) = (-1/4*u(:,6) + 4/3*u(:,5) - ...
      3*u(:,4) + 4*u(:,3) - 25/12*u(:,2))/dx;
  v_x(:,1) = (-1/4*v(:,5) + 4/3*v(:,4) - ...
      3*v(:,3) + 4*v(:,2) - 25/12*v(:,1))/dx;
  v_x(:,2) = (-1/4*v(:,6) + 4/3*v(:,5) - ...
      3*v(:,4) + 4*v(:,3) - 25/12*v(:,2))/dx;
  u_y(1,:) = (-1/4*u(5,:) + 4/3*u(4,:) - ...
      3*u(3,:) + 4*u(2,:) - 25/12*u(1,:))/dy;
  u_y(2,:) = (-1/4*u(6,:) + 4/3*u(5,:) - ...
      3*u(4,:) + 4*u(3,:) - 25/12*u(2,:))/dy;
  v_y(1,:) = (-1/4*v(5,:) + 4/3*v(4,:) - ...
      3*v(3,:) + 4*v(2,:) - 25/12*v(1,:))/dy;
  v_y(2,:) = (-1/4*v(6,:) + 4/3*v(5,:) - ...
      3*v(4,:) + 4*v(3,:) - 25/12*v(2,:))/dy;

  u_x(:,end) = -(-1/4*u(:,end-4) + 4/3*u(:,end-3) - ...
      3*u(:,end-2) + 4*u(:,end-1) - 25/12*u(:,end))/dx;
  u_x(:,end-1) = -(-1/4*u(:,end-5) + 4/3*u(:,end-4) - ...
      3*u(:,end-3) + 4*u(:,end-2) - 25/12*u(:,end-1))/dx;
  v_x(:,end) = -(-1/4*v(:,end-4) + 4/3*v(:,end-3) - ...
      3*v(:,end-2) + 4*v(:,end-1) - 25/12*v(:,end))/dx;
  v_x(:,end-1) = -(-1/4*v(:,end-5) + 4/3*v(:,end-4) - ...
      3*v(:,end-3) + 4*v(:,end-2) - 25/12*v(:,end-1))/dx;
  u_y(end,:) = -(-1/4*u(end-4,:) + 4/3*u(end-3,:) - ...
      3*u(end-2,:) + 4*u(end-1,:) - 25/12*u(end,:))/dy;
  u_y(end-1,:) = -(-1/4*u(end-5,:) + 4/3*u(end-4,:) - ...
      3*u(end-3,:) + 4*u(end-2,:) - 25/12*u(end-1,:))/dy;
  v_y(end,:) = -(-1/4*v(end-4,:) + 4/3*v(end-3,:) - ...
      3*v(end-2,:) + 4*v(end-1,:) - 25/12*v(end,:))/dy;
  v_y(end-1,:) = -(-1/4*v(end-5,:) + 4/3*v(end-4,:) - ...
      3*v(end-3,:) + 4*v(end-2,:) - 25/12*v(end-1,:))/dy;

  % Fix points that border a pore (ie. cut cells)
  % Allocate space for finding points whose stencil may cut through a
  % pore
  extPtsX = [];
  extPtsY = [];
  for k = 1:numel(radii)
    dist2 = (eX - centers(k,1)).^2 + (eY - centers(k,2)).^2;
    intPts = find(dist2 > (radii(k)-2*max(dx,dy))^2 & dist2 <= radii(k)^2);
    u_x(intPts) = 0;
    u_y(intPts) = 0;
    v_x(intPts) = 0;
    v_y(intPts) = 0;
    % Assign zero in the interior

    [newPtsY,newPtsX] = find(dist2 < ...
        (radii(k)+2*max(dx,dy))^2 & dist2 > radii(k)^2);
    extPtsX = [extPtsX; newPtsX]; 
    extPtsY = [extPtsY; newPtsY]; 
  end
%  figure(1)
%  plot(eX,eY,'rx')

  for k = 1:numel(extPtsX)
    j = extPtsY(k);
    i = extPtsX(k);
%    plot(eX(j,i),eY(j,i),'g.')
%    disp([eX(j,i) eY(j,i)])
%    pause
    % points whose second right neighbor to the right is inside a pore
    if any(((eX(j,i+2) - centers(:,1)).^2 + ...
          (eY(j,i) - centers(:,2)).^2) < radii.^2)
      u_x(j,i) = -(-1/4*u(j,i-4) + 4/3*u(j,i-3) - 3*u(j,i-2) + ...
          4*u(j,i-1) - 25/12*u(j,i))/dx;
      v_x(j,i) = -(-1/4*v(j,i-4) + 4/3*v(j,i-3) - 3*v(j,i-2) + ...
          4*v(j,i-1) - 25/12*v(j,i))/dx;
    end
    % points whose second left neighbor to the right is inside a pore
    if any(((eX(j,i-2) - centers(:,1)).^2 + ...
          (eY(j,i) - centers(:,2)).^2) < radii.^2)
      u_x(j,i) = (-1/4*u(j,i+4) + 4/3*u(j,i+3) - 3*u(j,i+2) + ...
          4*u(j,i+1) - 25/12*u(j,i))/dx;
      v_x(j,i) = (-1/4*v(j,i+4) + 4/3*v(j,i+3) - 3*v(j,i+2) + ...
          4*v(j,i+1) - 25/12*v(j,i))/dx;
    end
    % points whose second upper neighbor to the right is inside a pore
    if any(((eX(j,i) - centers(:,1)).^2 + ...
          (eY(j+2,i) - centers(:,2)).^2) < radii.^2)
      u_y(j,i) = -(-1/4*u(j-4,i) + 4/3*u(j-3,i) - 3*u(j-2,i) + ...
          4*u(j-1,i) - 25/12*u(j,i))/dy;
      v_y(j,i) = -(-1/4*v(j-4,i) + 4/3*v(j-3,i) - 3*v(j-2,i) + ...
          4*v(j-1,i) - 25/12*v(j,i))/dy;
    end
    % points whose second lower neighbor to the right is inside a pore
    if any(((eX(j,i) - centers(:,1)).^2 + ...
          (eY(j-2,i) - centers(:,2)).^2) < radii.^2)
      u_y(j,i) = (-1/4*u(j+4,i) + 4/3*u(j+3,i) - 3*u(j+2,i) + ...
          4*u(j+1,i) - 25/12*u(j,i))/dy;
      v_y(j,i) = (-1/4*v(j+4,i) + 4/3*v(j+3,i) - 3*v(j+2,i) + ...
          4*v(j+1,i) - 25/12*v(j,i))/dy;
    end
  end
end % order == 4
