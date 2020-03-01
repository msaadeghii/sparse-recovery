function [x,time] = lp_IRL1_Magic(A,b,p)

% Wriiten by Reza, based on the paper:
% ITERATIVELY REWEIGHTED ALGORITHMS FOR COMPRESSIVE SENSING
tic;
tol = 1e-3;
x = l1eq_pd(A' * b, A, [], b, tol);

d = inf;
Cntr = 0;
MaxItr = 50;
epsalg = 1;
while ( ( d > sqrt(epsalg) / 100 ) && (epsalg >= 2e-3) )

    y1 = x;
    
    Winv = diag( (abs(x)+ epsalg).^(1-p));
    
    temp = l1eq_pd(x, A * Winv, [], b, tol);
    x = (abs(x)+ epsalg).^(1-p) .* temp;
        
    d = norm(y1-x,2);% / norm(y1,2);

    Cntr = Cntr + 1;
    epsalg = epsalg / 10;
end;
time = toc;