%
% Ejercicio 4 pag 219
%
% Contrastamos Runge-Kutta de tres etapas
% y tablero de Butcher dado en el ejercicio
% con Euler explicito para una solucion donde 
% se obtiene orden 3.
%
% estimamos los ordenes de convergencia 
%  usando las rectas de pendiente
% Tablero de Butcher es 
%
% 0 |  0   0  0
% 1 |  1   0  0 
% 1 | 1/2 1/2 0
%---------------
%   | 3/6 1/6 2/6
%

t0=0; %Tiempo inicial
tf=1;  %Tiempo final
T=tf-t0;% Tiempo total
y0=1; %Dato inicial
%Datos para la solucion exacta
kk=-5;
% Solucion exacta en una particion fina
ttrue=t0:0.001:t0+T; %Particion fina

ytrue=y0*exp(kk*ttrue);
%
% Numero de calculos a realizar
%
M=15;
nP=zeros(1,M);% guarda el numero de puntos de la particion en cada calculo
errEuler=zeros(1,M); %guarda el error obtenido con Euler
errRK=zeros(1,M); %guarda el error obtenido con RK4
%
% Numero de puntos iniciales
%
N=100;
for j=1:M  
    nP(j)=N;% Se guarda el numero de puntos a usar 
    h=T/N; % Talla de la particion
    t=t0:h:t0+T; %Particion
    h=T/N;% talla de la particion
% Vector para Runge-Kutta orden 4
    yRK=zeros(1,N+1);% dimensionaliza t e y 
    yRK(1)=y0;
% Vector para Euler orden 1
    yEuler=zeros(1,N+1);% dimensionaliza t e y 
    yEuler(1)=y0;
% Solucion exacta en la particion 
    yt=y0*exp(kk*t);
for n=1:N
    %Calculo RK
    k1=mifun2(kk,t(n),yRK(n));
    k2=mifun2(kk,t(n)+h,yRK(n)+h*k1);
    k3=mifun2(kk,t(n)+h,yRK(n)+h/2*(k2+k1));
    yRK(n+1)=yRK(n)+h*(3*k1+1*k2+2*k3)/6;
    %Calculo Euler
    yEuler(n+1)=yEuler(n)+h*mifun2(kk,t(n),yEuler(n));
end 
figure(1);
plot(t,yEuler,'+-',t,yRK,'d-',ttrue,ytrue,'r-');
legend('Euler','RK','exacta','Location','Best');
title([' Modelo  con kk= ',num2str(kk),...
       ': RK vs Euler con N= ',num2str(N)]);
errEuler(j)=max(abs(yEuler-yt));
errRK(j)=max(abs(yRK-yt));
   pause(0.01); 
N=2*N; %duplicamos N
disp(['Muestra = ',num2str(j),' N= ',num2str(N),' Error Euler = ',num2str(errEuler(j)),...
    ' Error RK = ',num2str(errRK(j))])
pause(0.1);
end
%
% Visualizamos ahora los datos globales del calculo
%
figure(2)
plot(nP,log(errEuler),'-*',nP,log(errRK),'-+');
legend('LogErrEuler','LogErrRK','Location','Best');
title([' Modelo con kk= ',num2str(kk),...
       ': Decaimiento log errores Euler vs RK ']);
figure(3)
plot(log(nP),log(errEuler),'*',log(nP),-log(nP),'-',...
    log(nP),log(errRK),'bo',log(nP),-2*log(nP),'g-');
legend('Euler','-1','RK','-3','Location','Best');
title([' Modelo con kk= ',num2str(kk),...
       ': Ordenes convergencia Euler vs RK ']);
   