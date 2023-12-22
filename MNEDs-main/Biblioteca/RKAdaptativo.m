function [t,y] = RKAdaptativo(A, b, c, b_comparacion, f, y0, t0, T, N, tol, p)
  % Calcula la solucion de la EDO y'=f(t,y) por el metodo de Runge-Kutta

  % Parametros:
  % A tablero de Butcher para calcular cada pendiente
  % b vector con pesos para las pendientes (avance)
  % c vector suma fila i de A
  % b_comparacion vector para el error (mismo tamaño que b)
  % p orden del metodo de avance

  % f(t,y) funcion que define la EDO
  % (t0,y0) condicion inicial
  % [t0,t0+T] intervalo de la solucion
  % N número de puntos de la particion

  % Returns:
  % t vector con los valores de t en cada paso (la particion)
  % y vector con los valores de y en cada paso

    h = 0.5;
    t(1) = t0;
    Tfin = T + t0;
    y(1) = y0;
    n = 1;
    [m,m] = size(A);
    while((n <= N) && (t(n) < Tfin) && (h > eps))
      % RK para calcular pendientes
      q = 0; # pendiente promediada para avanzar
      k(1) = f(t(n), y(n)); # vector de pendientes
      q = q + k(1)*b(1);
      for j = 2 : m
        k(j) = f(t(n) + c(j)*h, y(n)+h*sum(A(j,1:j-1).*k(1:j-1)));
        q = q + k(j)*b(j);
      endfor
      dy = h*q;

      tau = abs(sum((b_comparacion-b).*k));
      % Aceptamos el paso dado y avanzamos
      if  (tau <= tol)
        y(n+1)=  y(n) + dy;
        t(n+1)=  t(n) + h;

        hNext=min(10*h,Tfin-t(n+1));
        n=n+1;
      % Fin de trabajo cuando se acepta el paso
      else
        %% hNext =h*s
        %% para ser usado si el error local de truncatura no es
        %% lo suficientemente pequeno
        %%
        hNext = 0.9*(tol/tau)^(1/p)*h;
        % Controlamos hNext
     end
    if (hNext< 0.1*h) % Que no sea menor que 0.1*h
      hNext = 0.1*h;
    end
    if (hNext > 10*h) % Que no sea mayor que 10*h
      hNext=10*h;
    end
    if (t(n) + hNext > Tfin)% Que no se sobrepase  t0+T
      hNext = Tfin - t(n);
      Nfin=n;
    end
    h = hNext;
    end

    if ((t(n) < Tfin) && (h <= eps))
      display("Error RK: no se pudo garantizar la tolerancia");
    endif

end

