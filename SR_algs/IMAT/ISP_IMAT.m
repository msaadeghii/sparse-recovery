function s=ISP_IMAT(A, x, maxiter, A_pinv)


% A_pinv=pinv(A);
s = A_pinv*x;
iter=1;

% IMATCS threshold parameters
alfa=0.1;T0=1000;
K=1:maxiter;
T=T0*exp(-alfa*(K-1));

% Main Loop

while iter <= maxiter %thr > tol
    
    thr=T(iter);
    
    s=s.*sign(max(abs(s)-thr,0));
    s = s - A_pinv*(A*s-x);
    
    iter=iter+1;
    
end