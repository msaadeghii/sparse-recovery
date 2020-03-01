function [x,t] = SCSA_FIT(A, b, eps1, eps2, c, lambda)

% This function implements the SCSA algorithm in the noisy case using the
% FIT instance of the algorithm.

% Inputs
% A: The sensing matrix
% b: The measurements vector
% eps1 : \epsilon_1 in the paper (outer loop stopping threshold)
% eps2 : \epsilon_2 in the paper (inner loop stopping threshold)
% c: The decay factor
% lambda: \lambda, the regularization parameter

% Outputs
% x: The recovered sparse vector
% t: Execution tim

tic;
[x,~] = BP_FISTA(A,b,lambda,eps2);

x_cur = x;
x_old = x_cur;

MaxIter = 5;
Cntr = 1;

sigma = 8 * c * max(abs(x));
A2 = A' * A;
Atb = A' * b;
maxeig = eigs(A2,1);
mu = 0.99 / ( 2 * maxeig  + lambda / sigma );

d1 = inf;

% while ( (Cntr <= MaxIter) && ( d1 >  eps1) )
while ( d1 >  eps1 )    
    y1 = x_cur;
    d2 = inf;
    
    
    Cntr2 = 0;
    y = x_cur;
    t_old = 1;
%     while (Cntr2 < 2000 && d2 >  eps2)
    while (d2 >  eps2)
        x_cur = shrink( y - 2 * mu * (A2 * y - Atb), lambda * sigma * mu , sigma);
        
        
        t_cur = (1 + sqrt(1 + 4 * t_old^2) ) / 2;
        y = x_cur + (t_old - 1) / t_cur * (x_cur - x_old);
        
        d2 = norm(x_cur - x_old, 2) / norm(x_old,2);
        
        x_old = x_cur;
        t_old = t_cur;
        Cntr2 = Cntr2 + 1;

    end;
    d1 = norm(y1 - x_cur,2) / norm(y1,2);
    sigma = sigma * c;
    mu = 0.99 / ( 2 * maxeig  + lambda / sigma );
%     Cntr = Cntr + 1;
    
    
end;
x = x_cur;
t = toc;