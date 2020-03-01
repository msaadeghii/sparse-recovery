function [x,t] = Chartrand_IT(A,b,lambda,eps,p)
tic;

A2 = A' * A;
Atb = A' * b;
maxeig = eigs(A2,1);
mu = 0.99 / ( 2 * maxeig );

[~,m] = size(A);
d = inf;
x = zeros(m,1);
Cntr = 0;
while (d > eps)
    y = x;
    x = shrink_p( x - 2 * mu * (A2 * x - Atb), lambda * mu,p);
    d = norm(y-x,2)/norm(y,2);
    Cntr = Cntr + 1;
end;
% Cntr
t = toc;