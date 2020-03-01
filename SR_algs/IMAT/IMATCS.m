% In the name of God
% IMAT_CS
%% general information
% This code implements the IMATCS compressed sensing recovery algorithm presented in
% the following papers:
% -M. Azghani and F.Marvasti,"iterative algorithms for random sampling and
% compressed sensing recovery", 2012.
% -M. Azghani and F. Marvasti, “Sparse signal processing,” in New Perspectives on Approximation and Sampling Theory.Springer,2014,pp.189–213.
%
% IMATCS algorithm solves the following problem:
%             min ||x||_0 
%             subject to y=Ax
% written by: Dr. Masoumeh Azghani (mazghani@sut.ac.ir)
% Version: 1.1
% Last modified: 9 August 2016.
% Affiliation: 
% Electrical Engineering Department, Sahand University of Technology
% Tabriz, Iran.
% For any problems, contact me at mazghani@sut.ac.ir
%% instructions of the code
% y is the measurement vector. 
% A is the measurement matrix.
% lambda is the relaxation parameter.
% T(.) is the thresholding operator which can be adjusted according to the application.
% An example for adjusting the thresholding operator for an image is given as:
% 
% itermax=50;alfa=0.1;T0=1000;
% K=1:itermax;
% T=T0*exp(-alfa*(K-1));

function [sols_IMATCS]=IMATCS(y,A,T)
itermax=length(T);
% itermax is the maximum number of iterations.
%% lambda calculation
D=A'*A;
EIG_A=eig(D);
maxEIG_A=max(EIG_A);
lambda=1.9/maxEIG_A;
%% renconstruction using IMAT_CS method
Y=A'*y;
X=zeros(size(Y));
for K=2:itermax
    X=lambda*Y+X-lambda*D*X;
    X=X.*(abs(X)>T(K));
end
sols_IMATCS=X;
