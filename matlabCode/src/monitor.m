classdef monitor
% Used for doing input/output of data.
% Can do plotting, write to files, compute error in area
% and length, and write to console

properties
verbose         % write data to console
saveData        % save data to the dat files and log file
axis            % axis of the plot
dataFile        % name of the file to write the data
logFile         % name of the file to write the log
Ninner          % number of points on inner boundaries
Nouter          % number of points on outer boundary
Ntracers        % number of tracers
nv              % number of inner boundaries
rtol            % relative tolerance for ode45
atol            % absolute tolerance for ode45

end % properties

methods

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function o = monitor(options,prams)
% monitor(X,Xwalls,options,prams) saves options and parameters
% needed by the class and also computes the initial error in
% area and length so that the errors in area and length
% can be monitored throughout the simulation.
% This is the constructor


o.verbose = options.verbose;
% write data to console
o.saveData = options.saveData;
% save messages to a log file
o.dataFile = options.dataFile;
% name of bin file for geometry
% Uses this file name to choose the file name for the tracers
o.logFile = options.logFile;
% name of log file
o.axis = options.axis;
% axis for the plotting
o.Ninner = prams.Ninner;
% number of points per inner boundary
o.Nouter = prams.Nouter;
% number of points per outer boundary
o.nv = prams.nv;
% number of inner boundaries


end % constructor: monitor

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function clearFiles(o)
% clearFiles() clears the log and bin file so that there is nothing from
% previous runs

fid = fopen(o.logFile,'w');
fclose(fid);
fid = fopen(o.dataFile,'w');
fclose(fid);
% delete the previous log and data files

end % clearFiles

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function welcomeMessage(o,options,prams)
% welcomeMessage(options,prams) writes specs from the simulation to the
% log file and console

o.writeStars
message = ['Num of points per inner boundary    ',...
    num2str(o.Ninner)];
o.writeMessage(message);
message = ['Num of points on the outer boundary ',...
    num2str(o.Nouter)];
o.writeMessage(message);
message = ['Num of exclusions                   ',...
    num2str(o.nv)];
o.writeMessage(message);
message = ['GMRES tolerance                     ',...
    num2str(prams.gmresTol,'%.0e')];
o.writeMessage(message);
if options.fmm
  message = 'FMM is being used';
else
  message = 'FMM is not being used';
end
o.writeMessage(message);

o.writeStars
o.writeMessage(' ');

end % welcomeMessage

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function writeStars(o)
% writeStars writes a message of stars to the console
% and the log file depending on verbose and saveData

messageStars = '*********************************************';
o.writeMessage(messageStars,'%s\n')

end % writeStars

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function writeMessage(o,message,format)
% function writeMessage(message,format) appends message 
% to o.fileName with format

if nargin == 2
  format = '%s\n';
end
% if user doesn't give format, take it to be string followed by a new
% line

if o.saveData
  fid = fopen(o.logFile,'a');
  fprintf(fid,format,message);
  fclose(fid);
end
% save to log file
if o.verbose
  disp(message)
end
% write to console if verbose==true


end % writeMessage


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function writeOutput(o,Xinner,Xouter,sigmaInner,sigmaOuter)
% writeOutput(Xinner,Xouter,sigmaInner,sigmaOuter) writes the geometry
% and the density function to .bin file.  This can be loaded later to do
% the tracer simulations without having to regenerate the density
% function

fid = fopen(o.dataFile,'w');
fwrite(fid,[o.Ninner;o.Nouter;o.nv],'double');
% write the number of points and boundaries
for k = 1:o.nv
  fwrite(fid,Xinner(1:o.Ninner,k),'double');
end
% write the x coordinate of the inner boundaries
for k = 1:o.nv
  fwrite(fid,Xinner(o.Ninner+1:end,k),'double');
end
% write the y coordinate of the inner boundaries
fwrite(fid,Xouter,'double');
% write the x and y coordinates of the outer boundary

for k = 1:o.nv
  fwrite(fid,sigmaInner(1:o.Ninner,k),'double');
end
% write the x coordinate of the inner densities
for k = 1:o.nv
  fwrite(fid,sigmaInner(o.Ninner+1:end,k),'double');
end
% write the y coordinate of the inner densities
fwrite(fid,sigmaOuter,'double');
% write the x and y coordinates of the outer density 

fclose(fid);
% close the bin file


end % writeOutput


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function writeTracers(o,eulerX,eulerY,u,v)
% writeTracers(eulerX,eulerY,u,v) writes the x and y coordinates of the
% Eulerian grid that is fixed and the x and y coordinates of the
% velocity at these points

[nx,ny] = size(eulerX);
% size of the output
fileName = [o.dataFile(1:end-8) 'Tracers.bin'];
% take the name of the Data file, strip the word Data, and add on the
% word Tracers

fid = fopen(fileName,'w');
fwrite(fid,[nx;ny],'double');
% write the sizes
for k = 1:ny
  fwrite(fid,eulerX(1:nx,k),'double');
end
% write the x coordinate at the fixed grid
for k = 1:ny
  fwrite(fid,eulerY(1:nx,k),'double');
end
for k = 1:ny
  fwrite(fid,u(1:nx,k),'double');
end
for k = 1:ny
  fwrite(fid,v(1:nx,k),'double');
end
fclose(fid);



end % writeTracers

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [Ninner,Nouter,nv,Xinner,Xouter,sigmaInner,sigmaOuter] = ...
    loadGeometry(o,fileName)

fid = fopen(fileName,'r');
val = fread(fid,'double');
fclose(fid);

Ninner = val(1);
Nouter = val(2);
nv = val(3);
val = val(4:end);

Xinner = zeros(2*Ninner,nv);
Xouter = zeros(2*Nouter,1);
sigmaInner = zeros(2*Ninner,nv);
sigmaOuter = zeros(2*Nouter,1);

istart = 1;
% start of a pointer to where everything is stored in val

for k = 1:nv
  iend = istart + Ninner - 1;
  Xinner(1:Ninner,k) = val(istart:iend);
  istart = iend + 1;
end
% x-coordinates of the inner boundaries

for k = 1:nv
  iend = istart + Ninner - 1;
  Xinner(Ninner+1:2*Ninner,k) = val(istart:iend);
  istart = iend + 1;
end
% y-coordinates of the inner boundaries

iend = istart + Nouter - 1;
Xouter(1:Nouter,1) = val(istart:iend);
istart = iend + 1;
% x-coordinates of the outer boundary

iend = istart + Nouter - 1;
Xouter(Nouter+1:2*Nouter,1) = val(istart:iend);
istart = iend + 1;
% y-coordinates of the outer boundary


for k = 1:nv
  iend = istart + Ninner - 1;
  sigmaInner(1:Ninner,k) = val(istart:iend);
  istart = iend + 1;
end
% x-coordinates of the inner densities

for k = 1:nv
  iend = istart + Ninner - 1;
  sigmaInner(Ninner+1:2*Ninner,k) = val(istart:iend);
  istart = iend + 1;
end
% y-coordinates of the inner densities 

iend = istart + Nouter - 1;
sigmaOuter(1:Nouter,1) = val(istart:iend);
istart = iend + 1;
% x-coordinates of the outer density

iend = istart + Nouter - 1;
sigmaOuter(Nouter+1:2*Nouter,1) = val(istart:iend);
istart = iend + 1;
% y-coordinates of the outer desnity


end % loadGeometry


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [nx,ny,eulerX,eulerY,u,v] = loadEuler(o,fileName)

fid = fopen(fileName,'r');
val = fread(fid,'double');
fclose(fid);

nx = val(1);
ny = val(2);
val = val(3:end);

eulerX = zeros(nx,ny);
eulerY = zeros(nx,ny);
u = zeros(nx,ny);
v = zeros(nx,ny);

istart = 1;
% start of a pointer to where everything is stored in val


for k = 1:ny
  iend = istart + nx - 1;
  eulerX(1:nx,k) = val(istart:iend);
  istart = iend + 1;
end

for k = 1:ny
  iend = istart + nx - 1;
  eulerY(1:nx,k) = val(istart:iend);
  istart = iend + 1;
end

for k = 1:ny
  iend = istart + nx - 1;
  u(1:nx,k) = val(istart:iend);
  istart = iend + 1;
end

for k = 1:ny
  iend = istart + nx - 1;
  v(1:nx,k) = val(istart:iend);
  istart = iend + 1;
end


end % loadEuler




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function plotData(o,Xinner,Xouter,eulerX,eulerY,u,v,xtra,ytra)

figure(1); clf; hold on;

quiver(eulerX,eulerY,u,v,'g')
plot(xtra,ytra,'r')
vec1 = [Xouter(1:end/2);Xouter(1)];
vec2 = [Xouter(end/2+1:end);Xouter(end/2)];
plot(vec1,vec2,'k')
fill(Xinner(1:end/2,:),Xinner(end/2+1:end,:),'k')
axis equal
axis(o.axis)



end % plotData


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function runMovie(o,Xinner,Xouter,xtra,ytra)

ntime = size(xtra,1);

figure(2);
for k = 1:ntime
  clf
  subplot(1,2,1)
  hold on
  vec1 = [Xouter(1:end/2);Xouter(1)];
  vec2 = [Xouter(end/2+1:end);Xouter(end/2)];
  plot(vec1,vec2,'k')
  fill(Xinner(1:end/2,:),Xinner(end/2+1:end,:),'k')

  plot(xtra(k,:),ytra(k,:),'r.')
  ax = o.axis;
  ax(4) = max(ytra(k,:))+2;
  ax(3) = ax(4) - 10;
  axis equal
  axis(ax)
  set(gca,'visible','off')

  subplot(1,2,2)
  hold on
  vec1 = [Xouter(1:end/2);Xouter(1)];
  vec2 = [Xouter(end/2+1:end);Xouter(end/2)];
  plot(vec1,vec2,'k')
  fill(Xinner(1:end/2,:),Xinner(end/2+1:end,:),'k')

  plot(xtra(k,:),ytra(k,:),'r.')
  ax = o.axis;
  ax(3) = min(ytra(k,:))-2;
  ax(4) = ax(3) + 10;
  axis equal
  axis(ax)
  set(gca,'visible','off')


  pause(0.01)
end





end % runMovie


end % methods


end % classdef

