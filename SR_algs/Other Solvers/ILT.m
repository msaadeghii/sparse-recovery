function [x,t] = ILT(A,b,lambda,eps)

% This function implements the ILT algorithm in D. Malioutov and A. Aravkin, 
% "Iterative log thresholding," in IEEE Int. Conf. Acoust. Speech Signal
% Process., 2014.

% Inputs
% A: The sensing matrix
% b: The measurements vector
% lambda: \lambda, the regularization parameter
% eps : stopping threshold

% Outputs
% x: The recovered sparse vector
% t: Execution tim

tic;

A2 = A' * A;
Atb = A' * b;
maxeig = eigs(A2,1);
mu = 0.99 / ( 2 * maxeig );


[~,m] = size(A);
d = inf;
x = zeros(m,1);
% Cntr = 0;
while (d > eps)
    y = x;
    x = shrink_log( x - 2 * mu * (A2 * x - Atb), lambda * mu);
    d = norm(y-x,2)/norm(y,2);
%     Cntr = Cntr + 1;
end;
% Cntr
t = toc;