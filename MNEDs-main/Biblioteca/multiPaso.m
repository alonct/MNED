function [t,y] = multiPaso(a, b, k, f, y0, t0, T, N)
  % Calcula la solución de la EDO y'=f(t,y) por el método de Runge-Kutta

  % Parámetros:
  % a coeficientes primer polinomio característico
    % a(j) coeficiente de y(n+j)
  % b coeficientes segundo polinomio característico
    % b(j) coeficiente de f_n+j

  % f(t,y) función que define la EDO
  % t0 tiempo inicial
  % y0 vector de k componentes con las condiciones iniciales:
  % y0(j+1) aproxima y(t0+j*h) para j = 0..k-1
  % [t0,t0+T] intervalo de la solución
  % N número de puntos de la partición

  % Returns:
  % t vector con los valores de t en cada paso (la partición)
  % y vector con los valores de y en cada paso

    h = T/N;
    t = t0:h:t0+T;
    y = zeros(1,N+1);
    y(1:k) = y0;
    for n = 1:N-k+1

        % Forma vectorial
        y(n+k) = - a*y(n:n+k-1)' + h*b*f(t(n:n+k-1),y(n:n+k-1))';

        % Forma secuencial
        %dy = 0; % incremento
        %for j = 0 : k-1
        %  dy = dy - a(j+1)*y(n+j) + h*b(j+1)*f(t(n+j),y(n+j));
        %endfor
        %y(n+k) = dy;
    end
end
