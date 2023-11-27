#function dydt=DBfun(tt,yy,a)
#dydt=a*(sin(tt)-yy);
#end

function dydt = DBfun(tt,yy,a)
    dydt = yy**2 - yy**3;
end