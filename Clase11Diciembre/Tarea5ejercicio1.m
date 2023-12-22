clear all;
%
% Difusion y transporte de dato discontinuo
%
format long;
%
% MATLAB indexa desde 1 los vectores
%
% Dimension de la particion M+1 xe(1),xe(2),...,xe(M+2) donde 
%   xe(1)=xIni, xe(M+2)=xFin y los nodos intermedios son 
%     xe(2),xe(3),...,xe(M+1)
%

M=200;
xIni = 0.0;
xFin =1.0;
h = (xFin - xIni)/(M+1);
xe=xIni:h:xFin;
h2=h*h;

% Uso de fuerza
%
f = zeros(1,M+2);
%
% Uso de fuerza inicial discontinua
% Muestra el uso de vectores
%
for j=1:M+2
    if ((xe(j)>0.3)&(xe(j)<0.5))
        f(j)=4;
    end
    if ((xe(j)>0.7)&(xe(j)<0.9))
        f(j)=1;
    end
end
%
% Coeficiente de difusion
%
k=1;

%
% u array solucion 
%
u1 = zeros(1,M+2); 


%    
% Construccion del sistema MxM
%    Diagonales de la matriz
%
dcent = 2*ones(M,1);  % d valores indexados de 1 a M
dinf = (-1)*ones(M-1,1);  % d-1 valores indexados de 1 a M-1
dsup = (-1)*ones(M-1,1);  % d-1 valores indexados de 1 a M-1
%
% Construimos matriz a partir de las diagonales principales
%
P=diag(dinf, -1) + diag(dcent, 0) + diag(dsup, 1);
for j=1:M
    P(j,j)=2+h2*(1+xe(j+1));
end
%
% termino independiente
%
d=zeros(M,1);


for j=1:M   % rango de 1 a M
    d(j)=h2*f(j+1);
end

z = P\d; 


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
plot(xe,u1,'.')
title([' Ejercicio 1 solucion con M = ',num2str(M),' puntos ']);

eig(P)


  
   

