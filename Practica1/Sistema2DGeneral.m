clear all;
T=5;
N=100;
x0=1;
y0=2;
h=T/N;
t=[0:h:T];
x0=4/3;
y0=2/3;
%
% Solucion exacta
%
xv=2*exp(-3*t)-exp(-39*t)+cos(t)/3;
yv=-exp(-3*t)+2*exp(-39*t)-cos(t)/3;
x=zeros(N+1,1);
y=zeros(N+1,1);
x(1)=x0;
y(1)=y0;
h=T/N;
ejes=[0 T -10 10];
%
% Matriz del sistema
%
A=[9 24;-24 -51];
%
% Matriz identidad
%
I=[1 0;0 1];
%
% Matriz de iteracion para Euler explicito
%
Me=I+h*A;
%
% Inicio computo con Euler explicito
%
v=[x(1);y(1)];
for i=1:N
    v=Me*v+h*[f1(t(i));f2(t(i))];
    x(i+1)=v(1);
    y(i+1)=v(2);
end
figure(1);
subplot(2,3,1),plot(t,x,'-o',t,xv,'-');
title([' Euler explicito, N= ',num2str(N)]);
legend('x','xv','location','Best');
%axis(ejes);
subplot(2,3,4),plot(t,y,'-o',t,yv,'-');
title([' Euler explicito, N= ',num2str(N)]);
legend('y','yv','location','Best');%axis(ejes);
%axis(ejes);
%
% Matriz de iteracion para Euler implicito
%
Mi=I-h*A;
%
% Inicio computo con Euler implicito
%
v=[x(1);y(1)];
for i=1:N  
    v=Mi\(v+h*[f1(t(i+1));f2(t(i+1))]);
    x(i+1)=v(1);
    y(i+1)=v(2);
end
subplot(2,3,2),plot(t,x,'*',t,xv,'-');
title([' Euler implicito, N= ',num2str(N)]);
legend('x','xv','location','Best');
axis(ejes);
subplot(2,3,5),plot(t,y,'*',t,yv,'-');
title([' Euler implicito, N= ',num2str(N)]);
legend('y','yv','location','Best');
axis(ejes);
%
% Matriz de iteracion para Crank-Nicolson
%
MCN=I-0.5*h*A;
MCNaux=I+0.5*h*A;
%
% Inicio computo con Crank-Nicolson
%
v=[x(1);y(1)];
for i=1:N   
    v=MCN\(MCNaux*v+0.5*h*[f1(t(i+1))+f1(t(i));f2(t(i+1))+f2(t(i))]);
    x(i+1)=v(1);
    y(i+1)=v(2);
end
subplot(2,3,3),plot(t,x,'+',t,xv,'-');
title([' Crank-Nicolson, N= ',num2str(N)]);
legend('x','xv','location','Best');
axis(ejes);
subplot(2,3,6),plot(t,y,'+',t,yv,'-');
title([' Crank-Nicolson, N= ',num2str(N)]);
legend('y','yv','location','Best');
axis(ejes);


