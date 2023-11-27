%
% estimamos el orden de convergencia con BDF3 
% usando las rectas de pendiente sobre el 
% problema y'=kk*cos(kk*t)*y que oscila periodicamente
% La solucion es 
%  y=y0*exp(sin(kk*t));
%
%Datos para la solucion exacta
t0=0; %Tiempo inicial
tf=20;  %Tiempo final
T=tf-t0;% Tiempo total
y0=2; %Dato inicial
%Datos para la solucion exacta
kk=7;

%
% Numero de calculos a realizar
%
M=13;
nP=zeros(1,M);% guarda el numero de puntos de la particion en cada calculo
errBDF3=zeros(1,M); %guarda el error obtenido con BDF3
%
% Numero de puntos iniciales
%
N=100;
for j=1:M  
    nP(j)=N;% Se guarda el numero de puntos a usar 
    h=T/N; % Talla de la particion
% Vectores para BDF3
    t=zeros(1,N+1);% dimensionaliza t
    y=zeros(1,N+1);% dimensionaliza y 
    yexact=zeros(1,N+1);% dimensionaliza yexact 
    %
    % Datos iniciales
    %
    t(1)=0;
    t(2)=h;
    t(3)=2*h;
    y(1)=y0*exp(sin(kk*t(1)));
    y(2)=y0*exp(sin(kk*t(2)));
    y(3)=y0*exp(sin(kk*t(3)));
    yexact(1)=y(1);
    yexact(2)=y(2);
    yexact(3)=y(3);
    
for n=3:N
    t(n+1)=n*h;
    y(n+1)=18/11*y(n)-9/11*y(n-1)+2/11*y(n-2);
    y(n+1)=y(n+1)/(1-6/11*h*kk*cos(kk*t(n+1)));
    yexact(n+1)=y0*exp(sin(kk*t(n+1)));
end 
figure(1);
plot(t,y,'x-',t,yexact,'r-');
legend('BDF3','exacta','Location','Best');
title([' BDF3 con N= ',num2str(N),' a= ',num2str(a)]);
%
errBDF3(j)=max(abs(yexact-y));
   pause(0.1); 
N=2*N; %duplicamos N
if (j>1)
   p= log(errBDF3(j-1)/errBDF3(j))/log(2);
disp(['Errores: N= ',num2str(N),' BDF3 = ',num2str(errBDF3(j)),' Orden ',num2str(p)])
end
pause(0.1);
end
%
% Visualizamos ahora los datos globales del calculo
%
figure(2)
plot(nP,log(errBDF3),'-*');
legend('LogErrBDF3','Location','Best');
title([' Modelo con a= ',num2str(a),...
       ': Decaimiento log errores BDF3 ']);
figure(3)
plot(log(nP),log(errBDF3),'*',log(nP),-3*log(nP),'-');
legend('BDF3','-3','Location','Best');
title([' Modelo con a= ',num2str(a),...
       ': Ordene de convergencia BDF3 ']);
   