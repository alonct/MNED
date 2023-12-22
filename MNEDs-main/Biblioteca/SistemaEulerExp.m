function [t, w] = SistemaEulerExp(A, F, d, w0, t0, T, N)
  % d dimensión

  h=T/N;
  t=[t0:h:T+t0];
  w=zeros(d,N);

  w(1:d,1) = w0;

  %
  % Matriz del sistema
  %

  %A

  %
  % Matriz identidad
  %

  I=eye(d);

  %
  % Matriz de iteracion para Euler explicito
  %

  Me=I+h*A;

  %
  % Inicio computo con Euler explicito
  %


  for i=1:N
      w(1:d,i+1)=Me*w(1:d,i)+h*[F(t(i),w(1:d,i))]';
  end

end
