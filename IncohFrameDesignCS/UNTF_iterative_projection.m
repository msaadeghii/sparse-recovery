function F=UNTF_iterative_projection(n, m, I)

% This function creats an (n*m) unit-norm tight-frame
% from a normal distribution using iterative projection
% onto the set of unit-norm frames and the set of
% tigh-frames. For more details on projection onto tight-frames, see [*].
%
% [*] J. A. Tropp, I. S. Dhillon, R. W. Heath Jr., and T. Strohmer, “Designing
%     structured tight frames via an alternating projection method,” 
%     IEEE Trans. Inf. Theory, vol. 51, no. 1, pp. 188–209, 2005.
%
% "I" denotes the number of iterations.

F = normc(randn(n,m));

for i=1:I
    % project onto tight-frames
    [U , ~, V] = svd(F);
    V=V(:,1:n);
    F=sqrt(m/n)*U*V';
    % project onto unit-norm frames
    F=normc(F);
end

return