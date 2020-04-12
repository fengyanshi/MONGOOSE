clear all
wave=load('waveform.out');
k=wave(:,1);
f=wave(:,2);
phi=wave(:,3);
nwave=length(wave);
xb=8.46;
tb=20.5;
aa=0.1084;
wc=5.53;
kc=3.2458;
x=[0:0.1:20];
t(1:201,1)=0:0.1:20;
dt=1.5;
[X T]=meshgrid(x,t);
[m n]=size(X);
ETA(1:m,1:n)=0;
for nw=1:nwave
eta(:,:,nw)=aa/nwave*cos(k(nw).*(X-xb)-2.*pi.*f(nw)*(T-tb));
ETA=ETA+eta(:,:,nw);
end
ETA=ETA*kc;
T=T*wc-dt*wc;
figure(1)
clf
subplot(511)
nx=69+1;
plot(T(:,nx),ETA(:,nx))
grid
xloc=x(nx);
title(['x=' num2str(xloc)])
subplot(512)
nx=54+1;
plot(T(:,nx),ETA(:,nx))
grid
xloc=x(nx);
title(['x=' num2str(xloc)])
subplot(513)
nx=38+1;
plot(T(:,nx),ETA(:,nx))
grid
xloc=x(nx);
title(['x=' num2str(xloc)])
subplot(514)
nx=23+1;
plot(T(:,nx),ETA(:,nx))
grid
xloc=x(nx);
title(['x=' num2str(xloc)])
subplot(515)
nx=7+1;
plot(T(:,nx),ETA(:,nx))
grid
xloc=x(nx);
title(['x=' num2str(xloc)])

break
figure(2)
clf
subplot(411)
nt=169;
plot(X(nt,:),ETA(nt,:))
grid
time=t(nt);
title(['t=' num2str(time)])
subplot(412)
nt=171;
plot(X(nt,:),ETA(nt,:))
grid
time=t(nt);
title(['t=' num2str(time)])
subplot(413)
nt=206;
plot(X(nt,:),ETA(nt,:))
grid
time=t(nt);
title(['t=' num2str(time)])
subplot(414)
nt=250;
plot(X(nt,:),ETA(nt,:))
grid
time=t(nt);
title(['t=' num2str(time)])




