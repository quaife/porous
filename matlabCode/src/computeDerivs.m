function [u_x,u_y,v_x,v_y] = computeDerivs(eX,eY,u,v);
% find gradient of velocity field using finite differences.  Want to use
% one-sided derivatives on the boundary and when the stencil is cut by a
% pore

order = 2;
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
end

if order == 2
  % Second-order at interior nodes
  u_x(:,2:end-1) = (u(:,3:end) - u(:,1:end-2))/2/dx;
  u_y(2:end-1,:) = (u(3:end,:) - u(1:end-2,:))/2/dy;
  v_x(:,2:end-1) = (v(:,3:end) - v(:,1:end-2))/2/dx;
  v_y(2:end-1,:) = (v(3:end,:) - v(1:end-2,:))/2/dy;
  % first-order finite difference for every point except the final
  % column/row

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
    disp(numel(radii) - k + 1)
    dist2 = (eX - centers(k,1)).^2 + (eY - centers(k,2)).^2;
    intPts = find(dist2 > 0.8*radii(k)^2 & dist2 <= radii(k)^2);
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

end

%clf
%plot(eX,eY,'r.')
%hold on;
%for k = 1:numel(extPtsX)
%  plot(eX(extPtsY(k),extPtsX(k)),eY(extPtsY(k),extPtsX(k)),'g.')
%end
%axis equal
%
