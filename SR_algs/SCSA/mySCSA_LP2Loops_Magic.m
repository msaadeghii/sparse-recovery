function [x,t,ItrNo] = mySCSA_LP2Loops_Magic(A, b, eps1, eps2, c,sigma_min)

% This function implements the SCSA algorithm in the no-noise case.
% l1_magic is used to solve the weighted l1 minimization problem.

% Inputs
% A: The sensing matrix
% b: The measurements vector
% eps1 : \epsilon_1 in the paper (outer loop stopping threshold)
% eps2 : \epsilon_2 in the paper (inner loop stopping threshold)
% c: The decay factor
% sigma_min: The minimum value for the paramter \sigma in the paper

% Outputs
% x: The recovered sparse vector
% t: Execution time
% ItrNo: Number of l1 minimization iterations done until reaching the final
% solution

tic;
opts.tol = 1e-3;
opts.print = 0;
x = yall1(A, b, opts);

MaxIter = 50;
Cntr = 1;

sigma = c * 8 * max(abs(x));

d1 = inf;

while ( ( d1 > eps1 ) && (Cntr <= MaxIter) && (sigma > sigma_min) )
    
    y1 = x;
    d2 = inf;
    
    while ( (d2 > eps2) && (Cntr <= MaxIter))
        y2 = x;
        
        w = (exp(-abs(x) / sigma));
        opts.weights=w;
        x = yall1(A, b, opts);

        
        d2 = norm(y2-x,2) / norm(y2,2);
        
        Cntr = Cntr + 1;
        
    end;
    
    d1 = norm(y1 - x,2) / norm(y1,2);
    
    sigma = sigma * c;
    
end;
% delta
ItrNo = Cntr;
t = toc;
