clear all;
close all;
addpath('./Biblioteca');

%% Test SistemaEulerExp
%% Sistema de trayectoria circular
%% x' = -y
%% y' = x
%% cond inicial (1,0)
%% solucion x=cost, y=sent

%% Tips: 1. Cuidado con las dimensiones de los parametros
%% y la salida de SistemaEulerExp
%% Asegurarse que:
%% - F devuelve un vector fila
%% - w0 es un vector columna
%% El valor devuelto w es una matriz d x N+1

test2 = false;

N0=5000;
t0 = 0;
T = 2*pi;

x0 = 1; %% Sistema de dim 2
y0 = 0;

d = 2;  % dimension
w0 = [x0 y0]; % dato inicial

A = [0 -1; 1 0]; % Matriz del sistema
% funcion del problema
F = @(t,w) [0 0]; % 0 porque el sistema es lineal

N = N0;
[t, ww] = SistemaEulerExp(A, F, d, w0', t0, T, N); %Metodo
w=ww'; % Trasponemos la salida para: filas = t, columnas = d
figure(1);
hold on;

for j = 1 : d
  plot(t, w(1:N+1,j), '-');
  k = 1;
endfor

legend('x', "y");
title(['Solucion modelo por Euler']);
% Se ve el coseno y el seno
hold off;

x = w(1:N+1,1);
y = w(1:N+1,2);

figure(2);
plot(x,y);
title("Trayectoria (x,y)");
% Una circunferencia

%% Test 2: solo para comprobar funcionamiento
if(test2)
  B = zeros(d); % Alternativa más fea
  G = @(t,w) [-w(2) w(1)];

  [t, ww] = SistemaEulerExp(B, G, d, w0', t0, T, N); %Metodo
  w=ww'; % Trasponemos la salida para: filas = t, columnas = d
  figure(3);
  hold on;

  for j = 1 : d
    plot(t, w(1:N+1,j), '-');
    k = 1;
  endfor

  legend('x', "y");
  title(['Solucion 2 modelo por Euler']);
  % Se ve el coseno y el seno
  hold off;

  x = w(1:N+1,1);
  y = w(1:N+1,2);

  figure(4);
  plot(x,y);
  title("Trayectoria 2 (x,y)");
  % Una circunferencia
endif

rmpath('./Biblioteca');
