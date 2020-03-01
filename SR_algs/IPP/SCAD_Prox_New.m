function y=SCAD_Prox_New(x,lam,a)

I1=abs(x) <= (2*lam);
I2=((2*lam) < abs(x)) & (abs(x) <= (a*lam));
I3=abs(x)> (a*lam);

y=sign(x).*max(abs(x)-lam,0).*I1+(((a-1)*x-sign(x)*(a*lam))/(a-2)).*I2+x.*I3;

end