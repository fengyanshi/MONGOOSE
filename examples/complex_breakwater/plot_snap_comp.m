
clear all;

fdir='./Results/';

data_status=load([fdir 'data_status.dat']);
imax=data_status(1,1); jmax=data_status(1,2); im1=imax-1;jm1=jmax-1;
prtdt=data_status(1,3);
ntype=data_status(1,4);

data_xi=load([fdir 'data_xi.dat']);
data_yj=load([fdir 'data_yj.dat']);


porous=load(['porous']);
obs=load('obs');
x=data_xi;
y=data_yj;

% draw domain
figure(1)
clf
porous(porous==1)=NaN;
pcolor(x,y,-porous),shading flat
hold on
obs(obs==1)=NaN;
pcolor(x,y,obs-1),shading flat
plot([0 400],[10 10],'k','LineWidth',2)
xlabel('x (m)')
ylabel('y (m)')
grid
print('-djpeg100','plots/complex_domain.jpg')

% draw grid
figure(2)
clf
[X,Y]=meshgrid(x,y);
skx=4;
sky=1;
line(X(1:sky:end,1:skx:end),Y(1:sky:end,1:skx:end))
line(X(1:sky:end,1:skx:end)',Y(1:sky:end,1:skx:end)')
xlabel('x (m)')
ylabel('y (m)')
print('-djpeg100','plots/complex_mesh.jpg')


figure(3)
clf
icount=0;
wid=6;
len=12;
set(gcf,'units','inches','paperunits','inches','papersize', [8 8],'position',[1 1 8 8],'paperposition',[0 0 8 8]);

nfile=[40 80 120];


for kt=1:length(nfile)
    file_num=sprintf('%.4d',nfile(kt));
    
data_f=load([fdir 'data_f.' file_num]);

subplot(length(nfile), 1, kt)
	 
 %data_f(obs<1)=NaN;
 pcolor(x,y,-data_f),shading flat
     hold on
     contour(x,y,-data_f,[-0.5 -0.5],'LineWidth',2,'Color','r','LineStyle','-')
title(['time = ' num2str(nfile(kt)*0.5) ' s '])
por=pcolor(x,y,porous),shading flat
     grid
     axis([0 400 0 12])
alpha(por,.5)
pcolor(x,y,obs),shading flat
end

print('-djpeg100','plots/complex_eta.jpg')

% k
figure(4)
clf
icount=0;
wid=6;
len=12;
set(gcf,'units','inches','paperunits','inches','papersize', [8 8],'position',[1 1 8 8],'paperposition',[0 0 8 8]);

nfile=[40 80 120];


for kt=1:length(nfile)
    file_num=sprintf('%.4d',nfile(kt));
    
data_f=load([fdir 'data_f.' file_num]);
data_k=load([fdir 'data_k.' file_num]);
data_u=load([fdir 'data_u.' file_num]);
data_v=load([fdir 'data_v.' file_num]);
a=size(data_k);  nmax=a(1,1)/jmax;

subplot(length(nfile), 1, kt)
	 
 pcolor(x,y,data_k),shading flat
caxis([0 0.1])
     hold on
     contour(x,y,-data_f,[-0.5 -0.5],'LineWidth',2,'Color','r','LineStyle','-')
title(['time = ' num2str(nfile(kt)*0.5) ' s '])
     grid
     axis([0 400 0 12])

% vector
skx=4;
sky=2;
scx=5;
scy=1;
quiver(X(1:sky:end,1:skx:end),Y(1:sky:end,1:skx:end),data_u(1:sky:end,1:skx:end)*scx,data_v(1:sky:end,1:skx:end)*scy,0)

por=pcolor(x,y,porous),shading flat
alpha(por,.5);
pcolor(x,y,obs+1),shading flat
end

print('-djpeg100','plots/complex_k.jpg')


% ----------
% k
figure(5)
clf
icount=0;
wid=6;
len=12;
set(gcf,'units','inches','paperunits','inches','papersize', [8 8],'position',[1 1 8 8],'paperposition',[0 0 8 8]);

nfile=[40 80 120];


for kt=1:length(nfile)
    file_num=sprintf('%.4d',nfile(kt));
    
data_f=load([fdir 'data_f.' file_num]);
data_k=load([fdir 'data_k.' file_num]);
data_u=load([fdir 'data_u.' file_num]);
data_v=load([fdir 'data_v.' file_num]);
data_p=load([fdir 'data_p.' file_num]);
a=size(data_k);  nmax=a(1,1)/jmax;

subplot(length(nfile), 1, kt)
	 
 pcolor(x,y,data_p),shading flat
 
%caxis([0 0.1])
     hold on
     contour(x,y,-data_f,[-0.5 -0.5],'LineWidth',2,'Color','r','LineStyle','-')
title(['time = ' num2str(nfile(kt)*0.5) ' s '])
     grid
     axis([0 400 0 12])

% vector
skx=4;
sky=2;
scx=5;
scy=1;
quiver(X(1:sky:end,1:skx:end),Y(1:sky:end,1:skx:end),data_u(1:sky:end,1:skx:end)*scx,data_v(1:sky:end,1:skx:end)*scy,0)

por=pcolor(x,y,porous),shading flat
alpha(por,.1);
pcolor(x,y,obs+20000),shading flat
end

print('-djpeg100','plots/complex_p.jpg')







 



 

