% Compute the numerical divergence of a velocity field with high-order
% stencils
dx = 0.05;
dy = 0.1;
[x,y] = meshgrid(-10:dx:10,-2:dy:2);

u = x./(x.^2 + (y-5).^2);
v = (y-5)./(x.^2 + (y-5).^2);
% some divergence free velocity field so we have something to test
% against
u_xTrue = (-x.^2 + (y-5).^2)./(x.^2 + (y-5).^2).^2;
v_yTrue = (+x.^2 - (y-5).^2)./(x.^2 + (y-5).^2).^2;


coeff_2 = [-1/2;0;1/2];
coeff_4 = [1/12;-2/3;0;2/3;-1/12];
coeff_6 = [-1/60;3/20;-3/4;0;3/4;-3/20;1/60];
coeff_8 = [1/280;-4/105;1/5;-4/5;0;4/5;-1/5;4/105;-1/280];

u_x = zeros(size(x));
v_y = zeros(size(y));
% space for the derivatives of the velocities

% START OF COMPUTING u_x
% 8th-order for points far enough in the interior
u_x(:,5:end-4) = (coeff_8(1)*u(:,1:end-8) + coeff_8(2)*u(:,2:end-7) + ...
                  coeff_8(3)*u(:,3:end-6) + coeff_8(4)*u(:,4:end-5) + ...
                  coeff_8(5)*u(:,5:end-4) + coeff_8(6)*u(:,6:end-3) + ...
                  coeff_8(7)*u(:,7:end-2) + coeff_8(8)*u(:,8:end-1) + ...
                  coeff_8(9)*u(:,9:end-0))/dx;

% 6th-order for next layer of points
ind = [4 size(u_x,2)-3];
u_x(:,ind) = (coeff_6(1)*u(:,ind-3) + coeff_6(2)*u(:,ind-2) + ...
              coeff_6(3)*u(:,ind-1) + coeff_6(4)*u(:,ind-0) + ...
              coeff_6(5)*u(:,ind+1) + coeff_6(6)*u(:,ind+2) + ...
              coeff_6(7)*u(:,ind+3))/dx;


% 4th-order for next layer of points
ind = [3 size(u_x,2)-2];
u_x(:,ind) = (coeff_4(1)*u(:,ind-2) + coeff_4(2)*u(:,ind-1) + ...
              coeff_4(3)*u(:,ind-0) + coeff_4(4)*u(:,ind+1) + ...
              coeff_4(5)*u(:,ind+2))/dx;

% 2th-order for next layer of points
ind = [2 size(u_x,2)-1];
u_x(:,ind) = (coeff_2(1)*u(:,ind-1) + coeff_2(2)*u(:,ind-0) + ...
              coeff_2(3)*u(:,ind+1))/dx;

% One-sided second-order difference for last layer of points
u_x(:,1) = (-1/2*u(:,3) + 2*u(:,2) - 3/2*u(:,1))/dx;
u_x(:,end) = (+1/2*u(:,end-2) - 2*u(:,end-1) + 3/2*u(:,end))/dx;
% END OF COMPUTING u_x


% START OF COMPUTING v_y
% 8th-order for points far enough in the interior
v_y(5:end-4,:) = (coeff_8(1)*v(1:end-8,:) + coeff_8(2)*v(2:end-7,:) + ...
                  coeff_8(3)*v(3:end-6,:) + coeff_8(4)*v(4:end-5,:) + ...
                  coeff_8(5)*v(5:end-4,:) + coeff_8(6)*v(6:end-3,:) + ...
                  coeff_8(7)*v(7:end-2,:) + coeff_8(8)*v(8:end-1,:) + ...
                  coeff_8(9)*v(9:end-0,:))/dy;

% 6th-order for next layer of points
ind = [4 size(v_y,1)-3];
v_y(ind,:) = (coeff_6(1)*v(ind-3,:) + coeff_6(2)*v(ind-2,:) + ...
              coeff_6(3)*v(ind-1,:) + coeff_6(4)*v(ind-0,:) + ...
              coeff_6(5)*v(ind+1,:) + coeff_6(6)*v(ind+2,:) + ...
              coeff_6(7)*v(ind+3,:))/dy;

% 4th-order for next layer of points
ind = [3 size(v_y,1)-2];
v_y(ind,:) = (coeff_4(1)*v(ind-2,:) + coeff_4(2)*v(ind-1,:) + ...
              coeff_4(3)*v(ind-0,:) + coeff_4(4)*v(ind+1,:) + ...
              coeff_4(5)*v(ind+2,:))/dy;

% 2th-order for next layer of points
ind = [2 size(v_y,1)-1];
v_y(ind,:) = (coeff_2(1)*v(ind-1,:) + coeff_2(2)*v(ind-0,:) + ...
              coeff_2(3)*v(ind+1,:))/dy;

% One-sided second-order difference for last layer of points
v_y(1,:) = (-1/2*v(3,:) + 2*v(2,:) - 3/2*v(1,:))/dy;
v_y(end,:) = (+1/2*v(end-2,:) - 2*v(end-1,:) + 3/2*v(end,:))/dy;
% END OF COMPUTING v_y


divergence = u_x + v_y;
surf(x,y,log10(abs(divergence)))
view(2); shading interp; axis equal;
colorbar
fprintf('Largest divergence in entire domain is   %6.2e\n',max(max(abs(divergence))))
fprintf('Largest divergence away from boundary is %6.2e\n',max(max(abs(divergence(5:end-4,5:end-4)))))



