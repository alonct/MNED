%% Ejercicio 2
%% @author Domingo Méndez
%% Ejecutar desde la carpeta contenedora
  % justo encima de Biblioteca

clear all;
close all;
addpath('./Biblioteca');

%% Inicializacion del problema

N0=10000;
t0 = 0;
T = 8;

x0 = -1; %% Sistema de dim 4
y0 = 0;
u0 = 0.1;
v0 = -0.1;

d = 4;  % dimension
w0 = [x0 y0 u0 v0]; % dato inicial
 % Coefs metodo
c = [0 0.5 0.5 1];
b = [1/6 1/3 1/3 1/6];
A = [0    0   0   0;
     0.5  0   0   0;
     0    0.5 0   0;
     0    0   1   0];
% funcion del problema
F = @(t,w) [w(3) w(4) -2*w(1)/(w(1)^2+w(2)^2) -2*w(2)/(w(1)^2+w(2)^2)];


N = N0;
[t, w] = SistemaRKExp(A, b, c, d, F, w0, t0, T, N); %Metodo

figure(1);
hold on;

for j = 1 : d
  plot(t, w(1:N+1,j), '-');
  k = 1;
endfor

legend('x', "y", "x'", "y'");
title(['Solucion modelo por RK-Clasico']);
hold off;

x = w(1:N+1,1);
y = w(1:N+1,2);
u = w(1:N+1,3); % u = x'
v = w(1:N+1,4); % v = y'

figure(2);
  plot(x,y);
  title("Ejercicio 2. Trayectoria (x,y) en [0,8]");
rmpath('./Biblioteca');

