clear all;
close all;
T=0.9;
a=-1;
b=3;
N=80;
h=(b-a)/N;
vel=1;
lambda=0.8;
k=vel*lambda*h;
x=a:h:b;
uold_ftbs=zeros(1,N+1);
unew_ftbs=uold_ftbs;


utrue=zeros(1,N+1);

for i=1:N+1
uold_ftbs(i)=u0(x(i));
end


tiempo=k;

for i=1:N+1
utrue(i)=uexact(tiempo,x(i),vel);
end
    %
    % Calulo FT-BS
    %
unew_ftbs(1)=bd(tiempo);
for i=2:N+1
 unew_ftbs(i)=uold_ftbs(i)-lambda*(uold_ftbs(i)-uold_ftbs(i-1));
 
end
   
plot(x,utrue,x,unew_ftbs,'*');
title(['Exacta vs FT-BS, M = ',num2str(N),', tiempo = ',...
    sprintf('%0.5g',tiempo),', lambda = ',sprintf('%0.5g',lambda)]);
legend('Exacta','FT-BS','Location','Best');
axis([-1 3 0 2]);

pause(1);
%
hold off;
uold_ftbs=unew_ftbs;


while (tiempo <= T)
    
tiempo=tiempo+k;
    %
    % Calulo FT-BS
    %
unew_ftbs(1)=bd(tiempo);
for i=2:N+1
 unew_ftbs(i)=uold_ftbs(i)-lambda*(uold_ftbs(i)-uold_ftbs(i-1));
end
  
for i=1:N+1
utrue(i)=uexact(tiempo,x(i),vel);
end
plot(x,utrue,x,unew_ftbs,'*');
title(['Exacta vs FT-BS, M = ',num2str(N),', tiempo = ',...
    sprintf('%0.5g',tiempo),', lambda = ',sprintf('%0.5g',lambda)]);
legend('Exacta','FT-BS','Location','Best');
axis([-1 3 0 2]);

    
hold off;
pause(1);
uold_ftbs=unew_ftbs;

end
     
