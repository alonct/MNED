clear all;
%
% Euler implicito y diferencias centradas en espacio
% Difusion de dato discontinuo
%
format long;
%
% MATLAB indexa desde 1 los vectores
%
% Dimension de la particion M+1 xe(1),xe(2),...,xe(M+2) donde 
%   xe(1)=xIni, xe(M+2)=xFin y los nodos intermedios son 
%     xe(2),xe(3),...,xe(M+1)
%

M=100;
xIni = 0.0;
xFin =1.0;
dx = (xFin - xIni)/(M+1);
%
% Talla de x es 1x(M+2)
%
xe =xIni:dx:xFin; % particion talla h, total de M+2 puntos

dx2=dx*dx;
%
% tiempo incial
%
tiempo=0.0;
%
% Valor inicial. 
% u0 array solucion en el paso de tiempo t^n
%
% Talla de u0 es 1x(M+2)
%
% Uso de dato inicial irregular
% Muestra el uso de vectores
%
u0 = zeros(1,M+2);

for j=1:M+2
    if ((xe(j)>0.3)&&(xe(j)<0.5))
        u0(j)=1;
    end
    if ((xe(j)>0.7)&&(xe(j)<0.9))
        u0(j)=4;
    end
end
%
% Coeficiente de difusion
%
k=0.1;
%
% u1 array solucion en el paso de tiempo t^{n+1}
%
u1 = zeros(1,M+2); 
%
% Tiempo final
%
T=1;
%
%
Nfin=100;
dt=T/Nfin;
mu=k*dt/dx2;
%    
% Construccion del sistema MxM
%    Diagonales de la matriz
%
dcent1 = (1+mu)*ones(M,1);  % d valores indexados de 1 a M
dinf1 = (-0.5*mu)*ones(M-1,1);  % d-1 valores indexados de 1 a M-1
dsup1 = (-0.5*mu)*ones(M-1,1);  % d-1 valores indexados de 1 a M-1
%
% Construimos matriz a partir de las diagonales principales
%
P=diag(dinf1, -1) + diag(dcent1, 0) + diag(dsup1, 1);
%    
% Construccion del sistema MxM
%    Diagonales de la matriz
%
dcent2 = (1-mu)*ones(M,1);  % d valores indexados de 1 a M
dinf2 = (0.5*mu)*ones(M-1,1);  % d-1 valores indexados de 1 a M-1
dsup2 = (0.5*mu)*ones(M-1,1);  % d-1 valores indexados de 1 a M-1
%
% Construimos matriz a partir de las diagonales principales
%
Q=diag(dinf2, -1) + diag(dcent2, 0) + diag(dsup2, 1);
%
% termino independiente
%
d=zeros(M,1);
%
% Dato inicial
%
plot(xe,u0,'+');
title(['Aprox con M = ',num2str(M),...
       ' Dato inicial ',num2str(tiempo)]);
axis([0 1 0 5]);
%pause;
%
% iteracion temporal
% 
for nt=1:Nfin
 %
 % Avanzamos en nivel de tiempo
 %
 tiempo=nt*dt;
%
%construccion termino independiente
%
% Asignacion de datos al termino 
% independendiente del sistema Ay=d
%
for j=1:M   % rango de 1 a M
    d(j)=u0(j+1);%+dt*f(j+1);
end
    %

z = P\(Q*d); 

for j=1:M
  u1(j+1)=z(j);
end
% Asignacion de los datos de contorno en x=xIni y en x=xFin
u1(1)=0.0; 
u1(M+2)=0.0; 
%
% Dibujamos la solucion aproximada
%
figure(1);
plot(xe,u1,'.-')
title([' kappa = ',num2str(k),...
       ', difusion con dx = ',num2str(dx),...
      ' dt ',num2str(dt),...
      ' tiempo ',num2str(tiempo)]);
  axis([0 1 0 5]);
  pause(.1);
%
%actualizamos
%
  for j=1:M+2
   u0(j)=u1(j);
  end
  pause;
end

  
   

