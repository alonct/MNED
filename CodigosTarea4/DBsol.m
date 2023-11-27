function ytt=DBsol(tt,a,y0)

aa2=a*a;
aux=1/(1+1/aa2); 
inva=1.0/a;   
ytt=y0*exp(-a*tt)+(sin(tt)-cos(tt)*inva+exp(-a*tt)*inva)*aux;
end