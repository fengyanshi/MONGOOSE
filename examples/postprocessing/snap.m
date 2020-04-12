
clear all;

fdir1='./Results/';
fdir10='./';
fdir2='../Cassion/Results/';
fdir20='../Cassion/';

data_status=load([fdir1 'data_status.dat']);
imax=data_status(1,1); jmax=data_status(1,2); im1=imax-1;jm1=jmax-1;
prtdt=data_status(1,3);
ntype=data_status(1,4);

data_xi=load([fdir1 'data_xi.dat']);
data_yj=load([fdir1 'data_yj.dat']);
obs1=load([fdir10 'obs']);
por1=load([fdir10 'porous']);
obs1(obs1==1)=NaN;
por1(por1==1)=NaN;

obs2=load([fdir20 'obs']);
por2=load([fdir20 'porous']);
obs2(obs2==1)=NaN;
por2(por2==1)=NaN;

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

data_f=load([fdir1 'data_f.' fnum]);
data_u=load([fdir1 'data_u.' fnum]);
data_v=load([fdir1 'data_v.' fnum]);
data_k=load([fdir1 'data_k.' fnum]);
data_p=load([fdir1 'data_p.', fnum]);
a=size(data_k);  nmax=a(1,1)/jmax;
isk=1;
colormap jet;

X=X1;
Y=Y1;
data_u(data_f<0.1)=NaN;
data_v(data_f<0.1)=NaN;
X(data_f<0.5)=NaN;
Y(data_f<0.5)=NaN;
data_f(obs1<1)=NaN;

clf
subplot(211)
pcolor(X(2:sk:jm1,2:sk:im1),Y(2:sk:jm1,2:sk:im1),data_p(2:sk:jm1,2:sk:im1)),shading interp
h=colorbar;
set(get(h,'ylabel'),'String','Pressure', 'rotation',270,'verticalalignment','bottom')
caxis([0 220000]);
axis([470 550 -3.3  26.6])
grid
hold on
%pcolor(X1,Y1,por1),shading interp
pcolor(X1,Y1,(obs1+1.5)),shading interp
time=sprintf('%5.2f',num*prtdt);
   title(['t=' time 's']);

ylabel('z (m)');
% -------------------------------------------
data_f=load([fdir2 'data_f.' fnum]);
data_u=load([fdir2 'data_u.' fnum]);
data_v=load([fdir2 'data_v.' fnum]);
data_k=load([fdir2 'data_k.' fnum]);
data_p=load([fdir2 'data_p.', fnum]);
a=size(data_k);  nmax=a(1,1)/jmax;
isk=1;
colormap jet;

X=X1;
Y=Y1;
data_u(data_f<0.1)=NaN;
data_v(data_f<0.1)=NaN;
X(data_f<0.5)=NaN;
Y(data_f<0.5)=NaN;
data_f(obs2<1)=NaN;

subplot(212)
pcolor(X(2:sk:jm1,2:sk:im1),Y(2:sk:jm1,2:sk:im1),data_p(2:sk:jm1,2:sk:im1)),shading interp
h=colorbar;
set(get(h,'ylabel'),'String','pressure', 'rotation',270,'verticalalignment','bottom')
caxis([0 220000]);
axis([470 550 -3.3 26.6])
grid
hold on
%pcolor(X1,Y1,por2),shading interp
pcolor(X1,Y1,(obs2+1.5)),shading interp

xlabel('x (m)');
ylabel('z (m)');


M(:,icount)=getframe(gcf);
pause(0.03)

end
break




 

