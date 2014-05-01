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


o.verbose = options.verbose;        % write data to console
o.saveData = true;      % save the data
o.axis = [0 5 -10 50];
o.dataFile = options.dataFile;
o.logFile = options.logFile;
o.Ninner = prams.Ninner;
o.Nouter = prams.Nouter;
o.Ntracers = prams.Ntracers;
o.nv = prams.nv;

if o.verbose
  fid = fopen(o.logFile,'w');
  fclose(fid);
  fid = fopen(o.dataFile,'w');
  fclose(fid);
end
% delete the previous log and data files



end % constructor: monitor

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function welcomeMessage(o,options,prams)

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
message = ['Num of tracers                      ',...
    num2str(o.Ntracers)];
o.writeMessage(message);
message = ['GMRES tolerance                     ',...
    num2str(prams.gmresTol,'%.0e')];
o.writeMessage(message);
message = ['ode45 relative tolerance            ',...
    num2str(prams.rtol,'%.0e')];
o.writeMessage(message);
message = ['ode45 absolute tolerance            ',...
    num2str(prams.atol,'%.0e')];
o.writeMessage(message);
message = ['ode45 time horizon                  ',...
    num2str(prams.T,'%2.1e')];
o.writeMessage(message);
message = ['no of time steps that ode45 returns ',...
    num2str(prams.ntime)];
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

fid = fopen(o.dataFile,'w');
fwrite(fid,[o.Ninner;o.Nouter;o.nv],'double');
for k = 1:o.nv
  fwrite(fid,Xinner(1:o.Ninner,k),'double');
end
for k = 1:o.nv
  fwrite(fid,Xinner(o.Ninner+1:end,k),'double');
end
fwrite(fid,Xouter,'double');
for k = 1:o.nv
  fwrite(fid,sigmaInner(1:o.Ninner,k),'double');
end
for k = 1:o.nv
  fwrite(fid,sigmaInner(o.Ninner+1:end,k),'double');
end
fwrite(fid,sigmaOuter,'double');
fclose(fid);




end % writeOutput

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function plotData(o)


end % plotData


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function writeTracers(o,z)
% writeData(X,sigma,ea,el,time,res) writes the position,
% tension, errors, and time to a binary file.  Matlab can later
% read this file to postprocess the data

oc = curve;
[x,y] = oc.getXY(z');

fid1 = fopen('tracersx.dat','w');
fid2 = fopen('tracersy.dat','w');
fprintf(fid1,'%4d\n',o.Ntracers);
fprintf(fid2,'%4d\n',o.Ntracers);
% write the number of points
for k = 1:size(x,1)
  fprintf(fid1,'%8.4e\n',x(k,:));
  fprintf(fid2,'%8.4e\n',y(k,:));
end
fclose(fid1);
fclose(fid2);


end % writeTracers
 


end % methods


end % classdef

