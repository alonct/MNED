
clear all;
%
% Uso del metodo BDF3 para la ecuacion de Dalhquist-Bjork
%
a=100;
t0=0; %Tiempo inicial
T=3;

N=100; % Numero de puntos maximo 
h=T/N;
% Vectores para BDF3
t=zeros(1,N+1);% dimensionaliza t
y=zeros(1,N+1);% dimensionaliza y 
yexact=zeros(1,N+1);% dimensionaliza yexact 



y0=1.0; %Dato inicial
%
% Datos iniciales
%
t(1)=0*h;
t(2)=1*h;
t(3)=2*h;
y(1)=DBsol(t(1),a,y0);
y(2)=DBsol(t(2),a,y0);
y(3)=DBsol(t(3),a,y0);
yexact(1)=y(1);
yexact(2)=y(2);
yexact(3)=y(3);


for n=3:N
    t(n+1)=n*h;
    y(n+1)=18/11*y(n)-9/11*y(n-1)+2/11*y(n-2)+6/11*h*a*sin(t(n+1));
    y(n+1)=y(n+1)/(1+6/11*h*a);
    yexact(n+1)=DBsol(t(n+1),a,y0);
end
figure(1);
plot(t,y,'x-',t,yexact,'r-');
legend('BDF3','exacta','Location','Best');
title([' BDF3 con N= ',num2str(N),' a= ',num2str(a)]);
%
% Calculo de los errores
%
errBDF3=abs(y-yexact);
figure(2);
semilogy(t,errBDF3,'o-');
legend('BDF3','Location','Best');
title([' Errores puntuales BDF3, N= ',num2str(N),' a= ',num2str(a)]);
