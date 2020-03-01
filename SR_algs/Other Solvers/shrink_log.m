function y = shrink_log(x,lambda)
eps = 1e-6;
x0 = (sqrt(2 * lambda) - eps) * ones(size(x));

y = 0.5 * (abs(x) > x0 ) .* ( x - eps * sign(x) + sign(x) .* sqrt( (x + eps * sign(x)).^2 - 2 * lambda )  ) ;