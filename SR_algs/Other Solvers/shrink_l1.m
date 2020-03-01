function z = shrink_l1(z0,lambda)

% z = zeros(size(z0));
zer = zeros(size(z0));

z = ( (abs(z0) - lambda) > zer) .* ( (abs(z0) - lambda) ) .* sign(z0);

% for i = 1 : length(z0)
%     if abs(z0(i)) >= lambda
%         z(i) = (abs(z0(i)) - lambda) * sign(z0(i));
%     else
%         z(i) = 0;
%     end;
% end;
