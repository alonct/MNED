%  Uso de script o guion
%
% Euler explicito 2D para sistema
%
% tiempo inicial
t0=0;

% tiempo final
T=20;
% Numero de puntos
N=1800;
% Paso
h=T/N;
% vector de tiempos
t=t0:h:t0+T;
hold off;
% vector de aproximaciones
x=zeros(1,N+1);
y=zeros(1,N+1);
% valor inicial
x(1)=0;
y(1)=2;

% calculo de Euler explicito
for s=1:N
    x(s+1)=x(s)+h*mifxBifurcacionHopf(x(s),y(s));
    y(s+1)=y(s)+h*mifyBifurcacionHopf(x(s),y(s));
end
% poner figura en primer plano
figure(1);
% dibuja los puntos calculados
subplot(2,1,1)
plot(t,x,'+',t,y,'-');
title(['curvas x e y. N= ',num2str(N),' T= ',num2str(T)]);
subplot(2,1,2)
plot(x,y,'-');
title(['Plano de fases. N= ',num2str(N),' T= ',num2str(T)]);
hold on;
%pause(2);
%end
