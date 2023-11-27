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

M=400;
xIni = 0.0;
xFin =1.0;
h = (xFin - xIni)/(M+1);
xe=xIni:h:xFin;
h2=h*h;

% Uso de fuerza
%
f = zeros(1,M+2);
%
for j=1:M+2
     f(j)=1;
end
%
% Coeficiente de difusion
%
k=0.05;
%
% Velocidad de transporte
%
av=-1.5;
%
% u array solucion 
%
u = zeros(1,M+2); 

%
% Numero de Peclet
%
Pe=av/k;
%    
% Construccion del sistema MxM
%    Diagonales de la matriz
%
dcent = 2*ones(M,1);  % d valores indexados de 1 a M
dinf = (-1-h*Pe/2)*ones(M-1,1);  % d-1 valores indexados de 1 a M-1
dsup = (-1+h*Pe/2)*ones(M-1,1);  % d-1 valores indexados de 1 a M-1
%
% Construimos matriz a partir de las diagonales principales
%
P=diag(dinf, -1) + diag(dcent, 0) + diag(dsup, 1);
%
% termino independiente
%
d=zeros(M,1);


for j=1:M   % rango de 1 a M
    d(j)=h2*f(j+1)/k;
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
title([' kappa = ',num2str(k),...
       ', vel. = ',num2str(av),...
       ', transporte y difusion con dx = ',num2str(h)]);



  
   

