function [x,t] = BP_Magic(A,b)
tic;
tol = 1e-3;
x = l1eq_pd(A' * b, A, [], b, tol);
t = toc;