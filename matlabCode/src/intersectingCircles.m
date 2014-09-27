function z = intersectingCirlces(N,rB,xB,yB,rW,xW,yW)

%N = 2^15;
omega = (0:N-1)*2*pi/N;

%rB = 3;
%xB = 0.2; yB = 0.2;
%rW = 3.2;
%xW = 1.8; yW = 1.2;

x0 = xW - xB;
y0 = yW - yB;
% move black circle to the origin

%figure(2); clf; 
%subplot(1,2,1);hold on
%fill(rB*cos(omega)+xB,rB*sin(omega)+yB,'k','edgecolor','w');
%fill(rW*cos(omega)+xW,rW*sin(omega)+yW,'w','edgecolor','w');
%axis equal

d = sqrt(x0^2 + y0^2);
% This is how far the white circle is from the black when everything is
% rotated to the x-axis

x = (d^2 - rW^2 + rB^2)/2/d;
y = sqrt((4*d^2*rB^2 - (d^2-rW^2+rB^2)^2))/2/d;


thetaOuter = zeros(2,1);
thetaOuter(1) = angle(x + 1i*y) + angle(x0 + 1i*y0);
thetaOuter(2) = angle(x - 1i*y) + angle(x0 + 1i*y0);

thetaInner = zeros(2,1);

thetaInner(1) = angle((rB*cos(thetaOuter(1))+xB - xW) + ...
      1i*(rB*sin(thetaOuter(1)) + yB - yW));
thetaInner(2) = angle((rB*cos(thetaOuter(2))+xB - xW) + ...
      1i*(rB*sin(thetaOuter(2)) + yB - yW));

%plot(rB*cos(thetaOuter(1)) + xB,rB*sin(thetaOuter(1)) + yB,'g.')
%plot(rB*cos(thetaOuter(2)) + xB,rB*sin(thetaOuter(2)) + yB,'r.')
%plot(rW*cos(thetaInner(1)) + xW,rW*sin(thetaInner(1)) + yW,'ko')
%plot(rW*cos(thetaInner(2)) + xW,rW*sin(thetaInner(2)) + yW,'bo') 

if thetaOuter(1) < 0
  thetaOuter(1) = thetaOuter(1) + 2*pi;
end
if thetaOuter(2) < 0
  thetaOuter(2) = thetaOuter(2) + 2*pi;
end
if thetaInner(1) < 0
  thetaInner(1) = thetaInner(1) + 2*pi;
end
if thetaInner(2) < 0
  thetaInner(2) = thetaInner(2) + 2*pi;
end

if thetaOuter(1) > thetaOuter(2)
  thetaOuter(1) = thetaOuter(1) - 2*pi;
end
arcOuter = rB*(thetaOuter(2) - thetaOuter(1));

if thetaInner(1) > thetaInner(2)
  thetaInner(1) = thetaInner(1) - 2*pi;
end
arcInner = rW*(thetaInner(2) - thetaInner(1));

Nouter = ceil(arcOuter/(arcOuter + arcInner)*N);
Ninner = N - Nouter;

theta1 = linspace(thetaOuter(1),thetaOuter(2),Nouter)';
theta2 = linspace(thetaInner(2),thetaInner(1),Ninner+2)';
theta2 = theta2(2:end-1);

z = [rB*cos(theta1)+xB ; rW*cos(theta2) + xW] + ...
    1i*[rB*sin(theta1)+yB ; rW*sin(theta2) + yW];

%subplot(1,2,2)
%plot(z,'b-o')
%hold on
%plot(z(1),'b.')
%axis equal
%pause

