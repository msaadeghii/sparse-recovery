function [x,time] = ReWL1_Magic(A, b, eps)

tol = 1e-3;
tic;
x = l1eq_pd(A' * b, A, [], b, tol);

d = inf;

while (d > eps)
    y = x;
    Winv = diag(abs(x) + eps);
    
    temp = l1eq_pd(x, A * Winv, [], b, tol);
    x = (abs(x) + eps) .* temp;
    d = norm(y - x,2) / norm(y,2);
end;
time = toc;