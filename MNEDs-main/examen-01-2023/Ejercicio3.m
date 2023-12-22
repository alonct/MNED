%% Ejercicio 3
%% @author Domingo Méndez
%% Ejecutar desde la carpeta contenedora
  % justo encima de Biblioteca

clear all;
close all;
addpath('./Biblioteca');

% Datos del problema
f = @(t,y) t.^3-2*t.*y;
g = @(t) 0.5*(t.^2-ones(1,length(t)))+exp(ones(1,length(t))-t.^2); %exacta
t0 = 1;
y0 = 1;

T = 4;

k1 = 2; % Número de pasos
k2 = 2;

% Multipaso
coef_a = 0; % Coefs metodo 1 alpha=0
b0 = 3-coef_a;
b1 = -(1+coef_a);
a_1 = [coef_a -(1+coef_a)];
b_1 = 0.5*[b0 b1];

coef_a = -5; % Coefs metodo 2 alpha = -5
a_2 = [coef_a -(1+coef_a)];
b_2 = 0.5*[b0 b1];

h = [0.1 0.05 0.025]; % varios h comparamos
for l = 1:length(h)
  clear 'tt0';
  N = round(T/h(l));
  tt0 = t0:h(l):t0+T;
  ytrue = g(tt0);
  y0_1 = ytrue(1:k1);
  y0_2 = ytrue(1:k2);

  % Resolvemos por los 2 métodos
  [t,y1] = multiPaso(a_1, b_1, k1, f, y0_1, t0, T, N);
  [t,y2] = multiPaso(a_2, b_2, k2, f, y0_2, t0, T, N);

  error1(l) = max(abs(ytrue-y1));
  error2(l) = max(abs(ytrue-y2));

  figure(1); % Se observa la divergencia del segundo metodo por dos motivos:
  % 1. El metodo no es 0-estable para alpha = -5
  % 2. El metodo es explicito pero el problema es mas o menos rigido:
  % solucion con una exponencial de negativos

  ttrue = t0:0.05:t0+T;
  yplot = g(ttrue);
  figure(1)
  plot(t, log(y1), '*b', t, log(y2), '*g', ttrue, log(yplot), '-r');
  xlabel('t');
  ylabel('log(y)');
  legend('alpha=0', 'alpha=-5', 'exacta');
  title(['Ejercicio 3, comparacion h=', num2str(h(l))]);

  figure(2)
  plot(t, log(y1), '*b', ttrue, log(yplot), '-r');
  xlabel('t');
  ylabel('log(y)');
  legend('alpha=0', 'exacta');
  title(['Ejercicio 3, a=0 h=', num2str(h(l))]);

  figure(3)
  plot(t, log(y2), '*g', ttrue, log(yplot), '-r');
  xlabel('t');
  ylabel('log(y)');
  legend('alpha=-5', 'exacta');
  title(['Ejercicio 3, a=-5 h=', num2str(h(l))]);

  disp('Pausa: espera 5s');
  pause(5);

  nP(l) = N;
endfor

pause(2)

%% Dibujamos grafica de los errores para comparar estabilidad
close all;
figure(4);
plot(log(nP), log(error1), '-b', log(nP), log(error2), '-g');
hold on;
legend('error a=0', 'error a=-5', 'Location','Best');
xlabel('log(nP)');
ylabel('log(err)');
title(['Ejercicio 3, comparacion errores a = 0 y -5 ']);
disp('En la figura observamos que es 0-estable para a=0 pero no para a=-5');
hold off;
rmpath('./Biblioteca');
