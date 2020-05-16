clear
close
%DESIRED WAVE FIELD HAS HM0=6M AND TP=14.4S
Tp=5.1;
gamma=3.3;
T=5100; %LENGTH OF THE RECORD FOR 500 WAVES
df=1/T;
g=9.81;
N=3001;
f=df:df:(N+1)*df;
w=2*pi*f;
Tper=2*pi./w;
fp=1/Tp; %PEAK FREQUENCY
wp=2*pi*fp;
alpha=0.09;%SHOULD BE MODIFIED TO OBTAIN DESIRED HMO
lf=length(f);
for i=1:lf
	if(f(i)>=fp)
		s=0.09;
	else
		s=0.07;
	end
	A(i)=exp(-1.25*(f(i)/fp)^-4);
	B(i)=exp(-(f(i)-fp)^2/(2*s^2*fp^2));
	C(i)=alpha*g^2/(w(i)^5);
	S(i)=C(i)*A(i)*gamma^(B(i));
	a(i)=sqrt(S(i)*2*df);
end
	plot(f,S)
	phase=-pi+2*pi*rand(1,lf);

dt=0.1;
srate=floor(1/dt);
t=0:dt:T;
H=0;
N=length(f);
for i=1:N
	H=H+a(i)*sin(w(i)*t+phase(i));
end
%figure
%plot(t,H)
%axis([0 t(length(t)) -8 8])
%ZERO_CROSSING COMPUTES HS, TM, TS FOR A GIVEN TIME SERIES. USED TO CHECK
%THAT THE GENERATED TIME SERIES DOES HAVE THE DESIRED HMO AND TP
[res1]=zero_crossing(H,srate) 
la=length(a);
out(1)=la;
out(2:la+1)=a;
out(la+2:2*la+1)=Tper;
out(2*la+2:3*la+1)=phase;
save 'spectral.txt' out -ascii

