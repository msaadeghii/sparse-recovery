function [x,t] = myBP_Magic(A,b)
tic;
opts.tol = 1e-3;
opts.print = 0;
x = yall1(A, b, opts);
t = toc;