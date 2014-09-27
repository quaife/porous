function z = appleGeom(N,rB,xB,yB,r1,x1,y1,r2,x2,y2)

omega = (0:N-1)*2*pi/N;

%rB = 3;
%xB = 0.2; yB = 0.2;
%rW = 3.2;
%xW = 1.8; yW = 1.2;


%figure(2); clf; 
%subplot(1,2,1);hold on
%fill(rB*cos(omega)+xB,rB*sin(omega)+yB,'k','edgecolor','w');
%fill(r1*cos(omega)+x1,r1*sin(omega)+y1,'w','edgecolor','r');
%fill(r2*cos(omega)+x2,r2*sin(omega)+y2,'w','edgecolor','w');
%axis equal


thetaOuter = zeros(4,1);
thetaInner = zeros(4,1);

x0 = x1 - xB;
y0 = y1 - yB;
% move black circle to the origin
d = sqrt(x0^2 + y0^2);
% This is how far the white circle is from the black when everything is
% rotated to the x-axis

x = (d^2 - r1^2 + rB^2)/2/d;
y = sqrt((4*d^2*rB^2 - (d^2-r1^2+rB^2)^2))/2/d;
thetaOuter(1) = angle(x + 1i*y) + angle(x0 + 1i*y0);
thetaOuter(2) = angle(x - 1i*y) + angle(x0 + 1i*y0);

x0 = x2 - xB;
y0 = y2 - yB;
% move black circle to the origin
d = sqrt(x0^2 + y0^2);
% This is how far the white circle is from the black when everything is
% rotated to the x-axis

x = (d^2 - r2^2 + rB^2)/2/d;
y = sqrt((4*d^2*rB^2 - (d^2-r2^2+rB^2)^2))/2/d;
thetaOuter(3) = angle(x + 1i*y) + angle(x0 + 1i*y0);
thetaOuter(4) = angle(x - 1i*y) + angle(x0 + 1i*y0);

thetaInner(1) = angle((rB*cos(thetaOuter(1))+xB - x1) + ...
      1i*(rB*sin(thetaOuter(1)) + yB - y1));
thetaInner(2) = angle((rB*cos(thetaOuter(2))+xB - x1) + ...
      1i*(rB*sin(thetaOuter(2)) + yB - y1));
thetaInner(3) = angle((rB*cos(thetaOuter(3))+xB - x2) + ...
      1i*(rB*sin(thetaOuter(3)) + yB - y2));
thetaInner(4) = angle((rB*cos(thetaOuter(4))+xB - x2) + ...
      1i*(rB*sin(thetaOuter(4)) + yB - y2));

%plot(x1,y1,'k.')
%plot(x2,y2,'k.')
plot(rB*cos(thetaOuter(1)) + xB,rB*sin(thetaOuter(1)) + yB,'g.')
plot(rB*cos(thetaOuter(2)) + xB,rB*sin(thetaOuter(2)) + yB,'r.')
plot(rB*cos(thetaOuter(3)) + xB,rB*sin(thetaOuter(3)) + yB,'m.')
plot(rB*cos(thetaOuter(4)) + xB,rB*sin(thetaOuter(4)) + yB,'b.')

if thetaOuter(1) < 0
  thetaOuter(1) = thetaOuter(1) + 2*pi;
end
if thetaOuter(2) < 0
  thetaOuter(2) = thetaOuter(2) + 2*pi;
end
if thetaOuter(3) < 0
  thetaOuter(3) = thetaOuter(3) + 2*pi;
end
if thetaOuter(4) < 0
  thetaOuter(4) = thetaOuter(4) + 2*pi;
end
if thetaInner(1) < 0
  thetaInner(1) = thetaInner(1) + 2*pi;
end
if thetaInner(2) < 0
  thetaInner(2) = thetaInner(2) + 2*pi;
end
if thetaInner(3) < 0
  thetaInner(3) = thetaInner(3) + 2*pi;
end
if thetaInner(4) < 0
  thetaInner(4) = thetaInner(4) + 2*pi;
end


if thetaOuter(1) > thetaOuter(4)
  thetaOuter(1) = thetaOuter(1) - 2*pi;
end
arcOuter1 = rB*(thetaOuter(4) - thetaOuter(1));

if thetaOuter(3) > thetaOuter(2)
  thetaOuter(3) = thetaOuter(3) - 2*pi;
end
arcOuter2 = rB*(thetaOuter(2) - thetaOuter(3));

if thetaInner(2) > thetaInner(1)
  thetaInner(2) = thetaInner(2) - 2*pi;
end
arcInner1 = r1*(thetaInner(1) - thetaInner(2));

if thetaInner(4) > thetaInner(3)
  thetaInner(4) = thetaInner(4) - 2*pi;
end
arcInner2 = r2*(thetaInner(3) - thetaInner(4));

totArc = arcOuter1 + arcOuter2 + arcInner1 + arcInner2;

Nouter1 = ceil(arcOuter1/totArc*N);
Nouter2 = ceil(arcOuter2/totArc*N);
Ninner1 = ceil(arcInner1/totArc*N);
Ninner2 = N - Nouter1 - Nouter2 - Ninner1;;

theta1 = linspace(thetaOuter(1),thetaOuter(4),Nouter1+1)';
theta1 = theta1(1:end-1);
theta2 = linspace(thetaInner(4),thetaInner(3)-2*pi,Ninner2+1)';
theta2 = theta2(1:end-1);
theta3 = linspace(thetaOuter(3),thetaOuter(2),Nouter2+1)';
theta3 = theta3(1:end-1);
theta4 = linspace(thetaInner(2),thetaInner(1)-2*pi,Ninner1+1)';
theta4 = theta4(1:end-1);

z = [rB*cos(theta1)+xB ; r2*cos(theta2)+x2; ...
     rB*cos(theta3)+xB ; r1*cos(theta4)+x1] + ...
 1i*[rB*sin(theta1)+yB ; r2*sin(theta2)+y2; ...
     rB*sin(theta3)+yB ; r1*sin(theta4)+y1];

%subplot(1,2,2)
%plot(z,'b-o')
%hold on
%plot(z(1),'b.')
%axis equal
%fill(real(z),imag(z),'y')

