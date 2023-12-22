function [t,y] = RKExp(A, b, c, f, y0, t0, T, N)
  % Calcula la solución de la EDO y'=f(t,y) por el método de Runge-Kutta

  % Parámetros:
  % A tablero de Butcher para calcular cada pendiente
  % b vector con pesos para las pendientes
  % c(i) = sum(A(i,:)) usualmente
    % vector con los avances en t para los muestreos de las pendientes
  % f(t,y) función que define la EDO
  % (t0,y0) condición inicial
  % [t0,t0+T] intervalo de la solución
  % N número de puntos de la partición

  % Returns:
  % t vector con los valores de t en cada paso (la partición)
  % y vector con los valores de y en cada paso

    h = T/N;
    t = t0:h:t0+T;
    y = zeros(1,N+1);
    y(1) = y0;
    [s,m] = size(A);
    k = zeros(1,m);
    for n = 1:N
        q = 0; # pendiente promediada para avanzar
        k(1) = f(t(n), y(n));
        q = q + k(1)*b(1);
        for i = 2 : m
          k(i) = f(t(n) + c(i)*h, y(n)+h*sum(A(i,1:i-1).*k(1:i-1)));
          q = q + k(i)*b(i);
        endfor
        y(n+1) = y(n) + h*q;
    end

end

