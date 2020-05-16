
clear all;

fdir='./Results/';

% input file numbers
nfile=[1:225];

% basic info
data_status=load([fdir 'data_status.dat']);
imax=data_status(1,1); jmax=data_status(1,2); im1=imax-1;jm1=jmax-1;
prtdt=data_status(1,3);
ntype=data_status(1,4);

data_xi=load([fdir 'data_xi.dat']);
data_yj=load([fdir 'data_yj.dat']);
obs=load('obs');
x=data_xi;
y=data_yj;

len=4;
wid=12;

set(gcf,'units','inches','paperunits','inches','papersize', [wid len],'position',[1 1 wid len],'paperposition',[0 0 wid len]);
% Set up file and options for creating the movie
vidObj = VideoWriter('movie.avi');  % Set filename to write video file
vidObj.FrameRate=10;  % Define the playback framerate [frames/sec]
open(vidObj);

for kt=1:length(nfile)

file_num=sprintf('%.4d',nfile(kt));
    
data_f=load([fdir 'data_f.' file_num]);
data_k=load([fdir 'data_k.' file_num]);

data_k(data_f<0.1)=NaN;

clf
pcolor(x,y,data_k),shading flat
hold on
caxis([0 0.1])
obs(obs==1)=NaN;
pcolor(x,y,obs-1),shading flat
     contour(x,y,-data_f,[-0.1 -0.1],'LineWidth',2,'Color','r','LineStyle','-')
ftim=sprintf('%.3f',nfile(kt)*0.25);
title(['time = ' ftim ' s '])
grid

axis([350 470 0 5])
pause(0.2)
    currframe=getframe(gcf);
    writeVideo(vidObj,currframe);  % Get each recorded frame and write it to filename defined above

end
close(vidObj)






 



 

