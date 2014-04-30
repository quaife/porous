load tracersx.dat;
load tracersy.dat;

Ntracers = tracersx(1);
tracersx = tracersx(2:end); tracersy = tracersy(2:end);
tracersx = reshape(tracersx,numel(tracersx)/Ntracers,Ntracers);
tracersy = reshape(tracersy,numel(tracersy)/Ntracers,Ntracers);

for k = 1:size(tracersx,1);
  clf; hold on;
  plot(tracersx(k,:),tracersy(k,:),'r.');
  plot(Xouter(1:end/2),Xouter(end/2+1:end),'k','linewidth',2);
  fill(Xinner(1:end/2,:),Xinner(end/2+1:end,:),'k');
  axis equal;

  ymax = max(tracersy(k,:))+3;
  ymin = ymax - 20;
  axis([-1 6 ymin ymax])
  set(gca,'xtick',[]);
  set(gca,'ytick',[]);
  set(gca,'visible','off');
  
%  fileName = ['image' num2str(k,'%04d') '.png'];
%  print('-dpng',fileName);

  pause()


end
