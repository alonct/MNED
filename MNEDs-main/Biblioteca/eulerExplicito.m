function [t,y, ztrue] = eulerExplicito(f, g, y0, t0, T, N)
  % Calcula la solución de la EDO y'=f(t,y) por el método de Euler Explícito
  % Y la dibuja

  % Parámetros:
  % f(t,y) función que define la EDO
  % g(t) la solución exacta de la ecuación
  % (t0,y0) condición inicial
  % [t0,t0+T] intervalo de la solución
  % N número de puntos de la partición

  % Returns:
  % t vector con los valores de t en cada paso (la partición)
  % y vector con los valores de y en cada paso
  % ztrue vector con la solución exacta en el intervalo
  % (en una partición fina h = 0.01 para dibujarla mejor)

    h = T/N;
    t = t0:h:t0+T;
    y = zeros(1,N+1);
    y(1) = y0;
    true = t0:0.01:t0+T;
    ztrue = g(true);
    for j=1:N
        y(j+1) = y(j) + h*f(t(j), y(j));
    end
    figure(1);
    plot(t,y,'*-',true,ztrue,'r-');

end

