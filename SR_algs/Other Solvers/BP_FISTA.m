function [x,t] = BP_FISTA(A,b,lambda,eps)

% This function implements FISTA algorithm for
% solving the unconstrained l1 minimization

% Inputs
% A: The sensing matrix
% b: The measurements vector
% lambda: \lambda, the regularization paramter
% eps: Stopping threshold

% Outputs
% x: The recovered sparse vector
% t: Execution time

tic;
A2 = A' * A;
Atb = A' * b;
maxeig = eigs(A2,1);
mu = 0.99 / ( 2 * maxeig );

[~,m] = size(A);
d = inf;
x_cur = zeros(m,1);
x_old = x_cur;
y = x_cur;
t_old = 1;
Cntr = 1;
while (d > eps)
    x_cur = shrink_l1( y - 2 * mu * (A2 * y - Atb), lambda * mu);
    t_cur = (1 + sqrt(1 + 4 * t_old^2) ) / 2;
    y = x_cur + (t_old - 1) / t_cur * (x_cur - x_old);
    
    d = norm(x_cur - x_old, 2) / norm(x_old,2);
    
    x_old = x_cur;    
    t_old = t_cur;
    Cntr = Cntr + 1;
end;
x = x_cur;
t = toc;