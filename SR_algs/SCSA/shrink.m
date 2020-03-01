function [ z ] = shrink(z0,alpha,sigma)

z = sign(z0) .* (sigma * lmbw0(-alpha * exp(-abs(z0)/sigma) / (sigma^2) ) + abs(z0));

y1 = 0.5*(z-z0).^2 + alpha *  (1 - exp(-abs(z)/sigma));
y2 = 0.5*z0.^2;

z = (y1 < y2) .* z;

end

