function [a,b,c]=abc(x1,y1,x2,y2,type)
% 1: /. .\ .|
% 2: ./ .\ .|
% 3: /. \.  |.
% 4. /. .\  |.
%top
if type ==1
if x1-x2~=0.
 s=(y2-y1)/(x2-x1);
 if s>=0
  a=-s;
  b=1.;
  c=s*x1-y1;
 else
  a=-s;
  b=1.;
  c=s*x1-y1;
 end
else
 a=1.;
 b=0.;
 c=-x1;
end
end

%left
if type ==2
if x1-x2~=0.
 s=(y2-y1)/(x2-x1);
 if s>=0
  a=s;
  b=-1.;
  c=-s*x1+y1;
 else
  a=-s;
  b=1.;
  c=s*x1-y1;
 end
else
 a=1.;
 b=0.;
 c=-x1;
end
end

%right
if type ==3
if x1-x2~=0.
 s=(y2-y1)/(x2-x1);
 if s>=0
  a=-s;
  b=1.;
  c=s*x1-y1;
 else
  a=s;
  b=-1.;
  c=-s*x1+y1;
 end
else
 a=-1.;
 b=0.;
 c=x1;
end
end

%bottom
if type ==4
if x1-x2~=0.
 s=(y2-y1)/(x2-x1);
 if s>=0
  a=-s;
  b=1.;
  c=s*x1-y1;
 else
  a=-s;
  b=1.;
  c=s*x1-y1;
 end
else
 a=-1.;
 b=0.;
 c=x1;
end
end
