
%
% El sistema de Lorenz se calcula mediante 
% Runge-Kutta de cuarto orden en el intervalo 
% [0,T] con h=T/N
%  con datos iniciales (x0,y0,z0)
%
% Se puede calcular la curva conforme se obtiene simplemente 
% dibujando las componentes de los vectores calculadas en cada 
% en cada paso de la interacion
 hold off;
 clear all;
 
 N=1000;
 T=10;
 x0=1;
 y0=2;
 z0=3;

h=T/N;
t=0:h:T;
x=zeros(1,N+1);
y=zeros(1,N+1);
z=zeros(1,N+1);
x(1)=x0;
y(1)=y0;
z(1)=z0;
%
% Aunque tengamos funciones de la forma f(x,y,z) y el argumento t 
% no sea necesario, es bueno ponerlo y escribir f(t,x,y,z) para 
% tener un codigo base que pueda servir para otros ejemplos.
for i=1:N
    k1x=fx(t(i),x(i),y(i),z(i));
    k1y=fy(t(i),x(i),y(i),z(i));
    k1z=fz(t(i),x(i),y(i),z(i));
    
    k2x=fx(t(i)+h/2,x(i)+h*k1x/2,y(i)+h*k1y/2,z(i)+h*k1z/2);
    k2y=fy(t(i)+h/2,x(i)+h*k1x/2,y(i)+h*k1y/2,z(i)+h*k1z/2);
    k2z=fz(t(i)+h/2,x(i)+h*k1x/2,y(i)+h*k1y/2,z(i)+h*k1z/2);
    
    k3x=fx(t(i)+h/2,x(i)+h*k2x/2,y(i)+h*k2y/2,z(i)+h*k2z/2);
    k3y=fy(t(i)+h/2,x(i)+h*k2x/2,y(i)+h*k2y/2,z(i)+h*k2z/2);
    k3z=fz(t(i)+h/2,x(i)+h*k2x/2,y(i)+h*k2y/2,z(i)+h*k2z/2);
    
    k4x=fx(t(i)+h,x(i)+h*k3x,y(i)+h*k3y,z(i)+h*k3z);
    k4y=fy(t(i)+h,x(i)+h*k3x,y(i)+h*k3y,z(i)+h*k3z);
    k4z=fz(t(i)+h,x(i)+h*k3x,y(i)+h*k3y,z(i)+h*k3z);
    
    x(i+1)=x(i)+h*(k1x+2*k2x+2*k3x+k4x)/6;
    y(i+1)=y(i)+h*(k1y+2*k2y+2*k3y+k4y)/6;
    z(i+1)=z(i)+h*(k1z+2*k2z+2*k3z+k4z)/6;
    
end
figure(1);
subplot(3,1,1),plot(t,x,t,y,t,z,'--');
title('Modelo de Lorenz')
legend('x','y','z');
hold off;
subplot(3,1,2),plot(x,y)
title('Efecto mariposa: Plano de Fases x-y')
hold off;
subplot(3,1,3),plot(x,z)
title('Efecto mariposa: Plano de Fases x-z')
hold off;


function valor=fx(tt,xx,yy,zz)
valor=-10*(xx-yy);
end

function valor=fy(tt,xx,yy,zz)
valor=-xx*zz+28*xx-yy;
end

function valor=fz(tt,xx,yy,zz)
valor=2.667*(xx*yy-zz);
end
