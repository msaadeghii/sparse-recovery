function [x,t] = BP_IST(A,b,lambda,eps)
% This function implements iterative soft thresholding algorithm for
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
x = zeros(m,1);

while (d > eps)
    y = x;
    x = shrink_l1( x - 2 * mu * (A2 * x - Atb), lambda * mu);
    d = norm(y-x,2)/norm(y,2);
end;
t = toc;