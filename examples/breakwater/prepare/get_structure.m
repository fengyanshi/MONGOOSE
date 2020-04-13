% program to generate structures
clear all;
inputname=['input_structure.txt'];
xy=load(inputname);
maxnum=floor(length(xy)/4);

outputname=['output_structure.txt'];
fn=outputname;
fid=fopen(fn,'wb');

clf
for num=1:maxnum
fnum=sprintf('%.2d',num);

ns=(num-1)*4+1;
ne=(num-1)*4+4;
xq=xy(ns:ne,1);
yq=xy(ns:ne,2);

[a(1),b(1),c(1)]= abc(xq(2),yq(2),xq(3),yq(3),1);

[a(2),b(2),c(2)]= abc(xq(1),yq(1),xq(2),yq(2),2);

[a(3),b(3),c(3)]= abc(xq(4),yq(4),xq(3),yq(3),3);

[a(4),b(4),c(4)]= abc(xq(1),yq(1),xq(4),yq(4),4);

qh(1)=1;
qh(2:4)=0;

% write out
line0=['3'];
dlmwrite(fn,line0,'-append')


for k=1:3
line1=[0.,a(k),0,b(k)];
line2=[0.,c(k),qh(k)];
dlmwrite(fn,line1,'-append')
dlmwrite(fn,line2,'-append')
end

% plot
Xq=[xq(1),xq(2),xq(3),xq(4),xq(1)];
Yq=[yq(1),yq(2),yq(3),yq(4),yq(1)];
fill(Xq,Yq,'y')
hold on
end
grid
xlabel('x(m)');
ylabel('y(m)');
print -djpeg100 structure.jpg



