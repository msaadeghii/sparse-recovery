function y=SCAD_Prox(X,lam,a)

n=length(X);
y=zeros(n,1);

for i=1:n
    x=X(i);
    if abs(x) <= (2*lam)
        y(i)=sign(x).*max(abs(x)-lam,0);
    elseif ((2*lam) < abs(x)) && (abs(x) <= (a*lam))
        y(i)=((a-1)*x-sign(x)*(a*lam))/(a-2);
    elseif abs(x)> (a*lam)
        y(i)=x;
    end
end