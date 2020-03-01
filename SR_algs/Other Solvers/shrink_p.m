function  z = shrink_p(s1,lambda,p)

z = zeros(size(s1));
ind = find( s1 ~= 0);
z(ind) = ( (abs(s1(ind)) - lambda * abs(s1(ind)).^(p - 1)) > 0) .* (abs(s1(ind)) - lambda * abs(s1(ind)).^(p - 1)) .* s1(ind) ./ abs(s1(ind));