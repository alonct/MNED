function [ uu] = uexact(tt,xx,vel)

if (abs(xx-vel*tt-1)<0.5)
uu=1;%cos(pi*(xx-tt))^2;
%uu=1-abs(xx-tt);
else
    uu=0;
end
end

