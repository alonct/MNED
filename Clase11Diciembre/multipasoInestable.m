
t0=0;
T=1;
N=20;
h=T/N
t=t0:h:t0+T;
true=t0:0.00001:t0+T;
ztrue=tan(true.^2);
y=zeros(N+1,1);
%
% Datos iniciales
%
y(1)=tan(t(1)^2);
y(2)=tan(t(2)^2);
y(3)=tan(t(3)^2);

for n=3:N
    y(n+1)=-1.5*y(n)+3*y(n-1)-0.5*y(n-2)+3*h*2*t(n)*(1+y(n)^2);
end
figure
plot(t,y,'*-',true,ztrue);
