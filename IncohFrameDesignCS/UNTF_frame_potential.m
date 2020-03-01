function F=UNTF_frame_potential(n, m, mu, I)

% This function creats an (n*m) unit-norm tight-frame
% from a normal distribution by minimizing the frame potential
% using the steepest descent algorithm. For more details, see [*].
%
% [*] J. J. Benedetto and M. Fickus, “Finite normalized tight frames,” 
%      Adv. Comput. Math., vol. 18, pp. 357–385, 2003.
%
% "mu" is the step-size parameter
% "I" denotes the number of iterations.

F=normc(randn(n,m));

for i=1:I
    F=F-mu*F*(F'*F);
    F=normc(F);
end

return