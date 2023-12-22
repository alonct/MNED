function [t,w] = SistemaRKExp(A, b, c, d, F, w0, t0, T, N)
  % Calcula la solución del sistema de EDOs y'=f(t,y)
  % por el método de Runge-Kutta explícito
  % con y de d dimensiones

  % Parámetros:
  % A tablero de Butcher para calcular cada pendiente
  % b vector con pesos para las pendientes
  % c dado por el método (sumas de las filas de A)
  % d dimensión del sistema

  % F(t,y) función que define la EDO (d-dimensional)
  % (t0,y0) condición inicial (y0 vector de tamaño d con el valor inicial)
  % T tamaño del intervalo: [t0,t0+T] intervalo de la solución
  % N número de puntos de la partición

  % Returns:
  % t vector con los valores de t en cada paso (la partición)
  % y matriz de tamaño [N+1 d] con la aproximación en cada paso


    h = T/N;
    t = t0:h:t0+T;
    y = zeros(1,N+1);

    w = zeros(N+1,d);
    w(1,:) = w0;
    [m,m] = size(A);
    for n = 1:N
        k = zeros(m,d);
        q = zeros(1,d); # pendiente promediada para avanzar
        k(1,:) = F(t(n), w(n,:));
        q += b(1)*k(1,:);
        y = zeros(1,d);
        for i = 2 : m
          y = w(n,:) + h*((A(i,1:i-1)*k(1:i-1,:)));
          k(i,:) = F(t(n) + c(i)*h, y);
          q += b(i)*k(i,:);
        endfor
        w(n+1,:) = w(n,:) + h*q;
    end
end
