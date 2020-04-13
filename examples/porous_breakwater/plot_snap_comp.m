
clear all;

fdir='./Results/';

data_status=load([fdir 'data_status.dat']);
imax=data_status(1,1); jmax=data_status(1,2); im1=imax-1;jm1=jmax-1;
prtdt=data_status(1,3);
ntype=data_status(1,4);

data_xi=load([fdir 'data_xi.dat']);
data_yj=load([fdir 'data_yj.dat']);


porous=load(['porous']);
x=data_xi;
y=data_yj;

% draw domain
figure(1)
clf
pcolor(x,y,porous),shading flat
hold on
plot([0 400],[10 10],'k','LineWidth',2)
xlabel('x (m)')
ylabel('y (m)')
print('-djpeg100','plots/porous_domain.jpg')

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
print('-djpeg100','plots/porous_mesh.jpg')


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
porous(porous==1)=NaN;
por=pcolor(x,y,porous),shading flat
     grid
     axis([0 400 0 12])
alpha(por,.5)
end

print('-djpeg100','plots/porous_eta.jpg')

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
	 
% data_k(obs<1)=NaN;
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

porous(porous==1)=NaN;
por=pcolor(x,y,porous),shading flat
alpha(por,.5);
end

print('-djpeg100','plots/porous_k.jpg')






 



 

