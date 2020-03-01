function mu=mucoh(F)
% computes the mutual coherence of a frame F
F=normc(F);   % first, normalize the frame
G=F'*F;  % compute the Gram matrix
A=G-eye(size(G));
mu=max(abs(A(:)));
return