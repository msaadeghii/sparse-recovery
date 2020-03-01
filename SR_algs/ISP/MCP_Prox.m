function y=MCP_Prox(x,lam,a)

I1=abs(x) <= (a*lam);
I2=abs(x) > (a*lam);

y=1/(1-1/a)*sign(x).*max(abs(x)-lam,0).*I1+x.*I2;

end