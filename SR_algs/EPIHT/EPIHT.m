function [xk,er,iter]=EPIHT(A, y, lam, tol, maxiter, beta, A_pinv, true_s)


% Initialization

x = A_pinv*y;
iter=1;
er=0;
xk=x;
xko=x;
diffs=inf;
L=eigs(A'*A,1);

sig=sqrt(2*lam/L);

% Main Loop
while  (iter < maxiter) %%(diffs > tol) &&
    xk1=xko;
    yk=xk+beta*(xk-xk1);
    if lam*sum(abs(yk)>1e-10)+1/2*norm(A*yk-y)^2 > lam*sum(abs(xk)>1e-10)+1/2*norm(A*xk-y)^2
       yk=xk;
    end
    grady=-A'*(y-A*yk);
    xki=yk-1/L*grady;
    xko=xk;
    xk=xki.*sign(max(abs(xki)-sig,0));
%     xk=mysoft(xki,sig^2);
%     diffs=norm(xko-xk)/(norm(xko));
    rel_res=norm(xk-true_s)/norm(true_s);
    er(iter)=rel_res;
    %     er(iter)=my_funcv(s,sigma);
    iter=iter+1;
end