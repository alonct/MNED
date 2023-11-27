clear all;
t0=0;
T=3;
y0=1;
N=20;
h=T/N;
a=10;
t=t0:h:t0+T;
y=zeros(N+1,1);
y(1)=y0;
true=t0:0.01:t0+T;
 %ztrue=true.^2/2;
 ztrue=exp(-a*true);
for j=1:N
    y(j+1)=y(j)+h*mifun(a,t(j),y(j));
end
figure(1);
plot(t,y,'*-',true,ztrue,'r-');

