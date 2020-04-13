
clear all;

%fdir='./Result_t8_a0.5_fine/';
fdir='./Results/';
fdir1='./Results/';

data_status=load([fdir 'data_status.dat']);
imax=data_status(1,1); jmax=data_status(1,2); im1=imax-1;jm1=jmax-1;
prtdt=data_status(1,3);
ntype=data_status(1,4);

data_xi=load([fdir 'data_xi.dat']);
data_yj=load([fdir 'data_yj.dat']);


obs=load([fdir 'data_ar.dat']);
%for ii=1:ntype
%por=load([fdir 'data_pr' num2str(ii) '.dat']);
%porous(:,:,ii)=por;
%end

%dummy=size(data_xi); im1=dummy(1);
%dummy=size(data_yj); jm1=dummy(1);

x=data_xi;
y=(data_yj-10)*2;


% objects
xobj1=[250 300 400 450 550 600 700 750];
xobj2=[250 300 400 450 550 600 700 750]*2;

%figure(1)
%clf;
%subplot 211

%-----------------------
%plot horizontal lines
%-----------------------
%y1=ones(size(x));

%x1=ones(size(y));


%nstart=input('input nstart=(530)')
%nend=input('input nend=(630)')

%nstart=1;
%nend=368;
icount=0;

set(gcf,'units','inches','paperunits','inches','papersize', [8 8],'position',[1 1 8 8],'paperposition',[0 0 8 8]);


sk=1;
st=5;
stv=3;
xst=11.;
xen=16;
nst=xst/(data_xi(2)-data_xi(1));
nen=xen/(data_xi(2)-data_xi(1));

%figure

nfile=[40 80 120 160];
nstt=40
%nfile=[nstt nstt+48 nstt+2*48 nstt+3*48];
clf

for kt=1:length(nfile)
    file_num=sprintf('%.4d',nfile(kt));
    
data_f=load([fdir 'data_f.' file_num]);
data_k=load([fdir 'data_k.' file_num]);
a=size(data_k);  nmax=a(1,1)/jmax;


subplot(4, 1, kt)
	 
for k=1:length(xobj2)
    data_f(:,xobj2(k)+1)=NaN;
        data_f(:,xobj2(k)+2)=NaN;
end

for k=1:length(xobj1)
    data_f1(:,xobj1(k)+1)=NaN;
        data_f1(:,xobj1(k)+2)=NaN;
end

	 x1=x(xobj2(1));
	 x2=x(xobj2(2));
	 y1=-1.1;
	 y2=1.1;
	 c = [0.6 0.8 1.0];
	 alpha=0.5;
	 h=fill([x1 x2 x2 x1 x1],[y1 y1 y2 y2 y1],c);
%	 set(h,'FaceColor',c,'FaceAlpha',alpha)
	hold on
	 x1=x(xobj2(3));
	 x2=x(xobj2(4));
	 c = [0.6 0.8 1.0];
	 alpha=0.5;
	 h=fill([x1 x2 x2 x1 x1],[y1 y1 y2 y2 y1],c);

	 x1=x(xobj2(5));
	 x2=x(xobj2(6));
     c = [0.6 0.8 1.0];
	 alpha=0.5;
	 h=fill([x1 x2 x2 x1 x1],[y1 y1 y2 y2 y1],c);

	 x1=x(xobj2(7));
	 x2=x(xobj2(8));

     c = [0.6 0.8 1.0];
	 alpha=0.5;
	 h=fill([x1 x2 x2 x1 x1],[y1 y1 y2 y2 y1],c);

     
	 if kt==4
	 xlabel('x (m)')
     end
	 ylabel('eta (m)')
	 set(gca,'LineWidth',1.2)
	 set(gcf,'renderer','opengl')

     hold on
     contour(x,y,data_f,[0.5 0.5],'LineWidth',2,'Color','r','LineStyle','--')
     contour(x11,y11,data_f1,[0.5 0.5],'LineWidth',2,'Color','b','LineStyle','-')

title(['time = ' num2str(nfile(kt)*0.5) ' s '])
     grid
     axis([0 800 -1 1])
end






 



 

