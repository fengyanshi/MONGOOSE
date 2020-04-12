
clear all;

fdir='./Results/';

data_status=load([fdir 'data_status.dat']);
imax=data_status(1,1); jmax=data_status(1,2); im1=imax-1;jm1=jmax-1;
prtdt=data_status(1,3);
ntype=data_status(1,4);

data_xi=load([fdir 'data_xi.dat']);
data_yj=load([fdir 'data_yj.dat']);
obs=load([fdir 'data_ar.dat']);
por=load('porous');
obs(obs==1)=NaN;
por(por==1)=NaN;

x=data_xi;
y=data_yj;

y1=ones(size(x));

x1=ones(size(y));

[X1,Y1]=meshgrid(x,y);
[mys mxs]=size(X1);

X=X1;
Y=Y1;
pcolor(X,Y,obs),shading flat
hold on
pcolor(X,Y,por),shading flat

xlabel('x (m)');
ylabel('z (m)');

