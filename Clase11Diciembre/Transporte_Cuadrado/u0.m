function [uu] = u0(xx)
if (abs(xx-1)<0.5)
uu=1;%cos(pi*xx)^2;
%uu=1-abs(xx);
else
    uu=0;
end

