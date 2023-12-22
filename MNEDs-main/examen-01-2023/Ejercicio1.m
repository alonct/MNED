%% Ejercicio 1
%% @author Domingo Méndez
%% Ejecutar desde la carpeta contenedora
  % justo encima de Biblioteca

clear all;
close all;
addpath('./Biblioteca');

%% Orden del método de Runge-Kutta

c = [0 0.5 0.75];
b = (1/9)*[2 3 4];

A = [0 0 0;
0.5 0 0;
0 0.75 0];

% Datos del problema
f = @(t,y) t.^3-2*t.*y;
y0 = 1;

% Sol exacta
g = @(t) 0.5*(t.^2-ones(1,length(t)))+exp(ones(1,length(t))-t.^2);

T = 4;
t0 = 1;

N = 10;
M = 16;

for k = 1 : M % Iteramos para calcular errores
  nP(k) = N;
  [t,y] = RKExp(A, b, c, f, y0, t0, T, N); % Metodo
  ytrue = g(t);
  error(k) = max(abs(ytrue-y));
  figure(1);
  ttrue = t0:0.05:t0+T;
  yplot = g(ttrue);
  plot(t, y, '*b', ttrue, yplot, '-r');
  title(['Solucion RK para Npuntos =', num2str(N)]);
  pause(0.2);

  N = round(1.5*N);
endfor

  figure(2);  % Mostramos la recta pendiente del orden y los errores
  plot(log(nP), log(error), '-b');
  hold on;
  j = M;
  p = 3;
  plot(log(nP(1:j)), -p*log(nP(1:j))+p*log(nP(j))+log(error(j)), '-r');
  legend('errorRK', ['recta p=' num2str(p)], 'Location','Best');
  xlabel('log(nP)');
  ylabel('log(err)');
  title(['Ejercicio 1, orden ' num2str(p) ]);

rmpath('./Biblioteca');
