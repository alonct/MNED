function dydt=mifun2(kk,tt,yy)
%dydt=kk*cos(kk*tt)*yy;
dydt=kk*yy;
%dydt=5*exp(5*tt)*(yy-tt)^2+1;
%dydt=-200*tt*yy^2;
end