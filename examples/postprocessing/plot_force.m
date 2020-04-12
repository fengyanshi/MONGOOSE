clear all
f=load('force.dat');
t1=f(:,1);
f1=f(:,2);
b1=f(:,3);
u1=f(:,4);
d1=f(:,5);
rl1=f(:,6);
rr1=f(:,7);

f=load('../Cassion/force.dat');
t2=f(:,1);
f2=f(:,2);
b2=f(:,3);
u2=f(:,4);
d2=f(:,5);
rl2=f(:,6);
rr2=f(:,7);

figure

plot(t1,b1,'b-',t2,b2,'r--','LineWidth',2)
legend('parapet','caission')
axis([700 890 800 2300])
grid
xlabel('t (s)')
ylabel('F ')
title('seaward force')

figure
plot(t1,rr1,'b-',t2,rr2,'r--','LineWidth',2)
legend('parapet','caission')
axis([700 890 10000 35000])
grid
xlabel('t (s)')
ylabel('M ')
title('landward moment')
break

figure
plot(t1,f1,'b-',t2,f2,'r--','LineWidth',2)
legend('parapet','caission')
axis([700 890 800 2300])
grid
xlabel('t (s)')
ylabel('F ')
title('forward force')

figure
plot(t1,rl1,'b-',t2,rl2,'r--','LineWidth',2)
legend('parapet','caission')
axis([700 890 -30000 -5000])
grid
xlabel('t (s)')
ylabel('rotional M ')
title('rotional moment')

figure
plot(t1,(u1-d1),'b-',t2,(u1-d1),'r--','LineWidth',2)
legend('parapet','caission')
axis([700 890 1950 2500])
grid
xlabel('t (s)')
ylabel('lift force ')
title('lift force')

