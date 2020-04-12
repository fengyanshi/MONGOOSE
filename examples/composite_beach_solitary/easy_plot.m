
clear all;

fdir='../Result11/';

data_status=load([fdir 'data_status.dat']);
imax=data_status(1,1); jmax=data_status(1,2); im1=imax-1;jm1=jmax-1;
prtdt=data_status(1,3);
ntype=data_status(1,4);

data_xi=load([fdir 'data_xi.dat']);
data_yj=load([fdir 'data_yj.dat']);
obs=load([fdir 'data_ar.dat']);
%por=load('porous');
obs(obs==1)=NaN;
%por(por==1)=NaN;

x=data_xi;
y=data_yj;

y1=ones(size(x));

x1=ones(size(y));

[X1,Y1]=meshgrid(x,y);
[mys mxs]=size(X1);

X=X1;
Y=Y1;

nstart=input('input nstart=(530)')
nend=input('input nend=(630)')

icount=0;

sk=1;
st=5;
stv=3;
xst=11.;
xen=16;
nst=xst/(data_xi(2)-data_xi(1));
nen=xen/(data_xi(2)-data_xi(1));

%figure

for num=nstart:1:nend
icount=icount+1;
fnum=sprintf('%.4d',num);

data_f=load([fdir 'data_f.' fnum]);
data_u=load([fdir 'data_u.' fnum]);
data_v=load([fdir 'data_v.' fnum]);
data_k=load([fdir 'data_k.' fnum]);

a=size(data_k);  nmax=a(1,1)/jmax;

isk=1;
 
%colormap cool;
colormap jet;


X=X1;
Y=Y1;

data_u(data_f<0.5)=NaN;
data_v(data_f<0.5)=NaN;
X(data_f<0.1)=NaN;
Y(data_f<0.1)=NaN;
data_f(obs<1)=NaN;
clf
colordef white
%axes('position',[.05 .1 .9 .35])
%subplot(211)
pcolor(X(2:sk:jm1,2:sk:im1),Y(2:sk:jm1,2:sk:im1),data_k(2:sk:jm1,2:sk:im1)),shading interp
%caxis([0 0.01])
caxis([0. 0.3])
axis([20 23.2 0.05 0.51])
grid
hold on
%pcolor(X,Y,por),shading interp
%pcolor(X,Y,obs),shading interp
%v=[0.02:0.04:1.5];
%hd=contour(X,Y,data_k,v,'y');

xf=[15.0,19.4,22.3,23.2,23.2,15.0];
yf=[0,0.083, 0.1023, 0.1715,0,0.0];

%text(14.8, 0.6,'bubble size: 5E-4')
xlabel('x (m)');
ylabel('z (m)');
facx=0.1;
facy=0.1;
quiver(X(2:stv:jm1,2:stv:im1),Y(2:stv:jm1,2:stv:im1),data_u(2:stv:jm1,2:stv:im1)*facx,data_v(2:stv:jm1,2:stv:im1)*facy,0,'Color','y')
fill(xf,yf,[.53 .44 0.08],'EdgeColor','None')

  time=sprintf('%5.2f',num*prtdt);
   title(['t=' time 's']);

xlabel('x (m)');
ylabel('z (m)');
M(:,icount)=getframe(gcf);
pause(0.03)

end
break




 

