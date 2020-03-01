function [F,mu]=IFD_AMPM(F0, muF, No, opts)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% This function implements the incoherent frame desing algorithm, dubbed IFD-AMPM,
% proposed in the following paper:
%
% [*] M. Sadeghi and M. Babaie-Zadeh, 
%     “Incoherent Unit-norm Frame Design via an Alternating Minimization Penalty Method,” 
%     IEEE Signal Processing Letters, 2016.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  Syntax:
%
%               [F,mu] = IFD_AMPM(F0, muF, No);                         defualt version
%               [F,mu] = IFD_AMPM(F0, muF, No, opts);               full version
%
% %%%%%%%%%%%%%%%%%%%% INPUTS %%%%%%%%%%%%%%%%%%%%%%%
%
%  Required:
%
%                       F0:                                  initial frame (can be set using "UNTF_iterative_projection") 
%                       muF:                              step-size of the F-update problem
%                       No:                                 number of iterations
% Optional:
%
%                        opts.c:                            threshold decaying factor (default: c=0.9)
%                        opts.alpha0:                 factor of the initial value for alpha (default: alpha0=500)
%                        opts.wF :                       weighting constant of F-update (default: wF=0.85)
%                        opts.wQ :                      weighting constant of Q-update (default: wQ=0.85)
%                        opts.Ni :                        number of inner-loop iterations (default: Ni=15)
%                        opts.T :                          number of F-update iterations (default: T=3)
%                                                                                                                          
% %%%%%%%%%%%%%%%%%%%% OUTPUTs %%%%%%%%%%%%%%%%%%%%%
%
% The algorithm returns "F", an incoherent frame, together with "mu", a vector
% containing the mutual coherence values along the iterations of the algorithm.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%     Mostafa Sadeghi                                                                                                   %%%
%%%     Electrical Engineering Department,                                                                 %%%
%%%     Sharif University of Technology,                                                                       %%%
%%%     Tehran, IRAN                                                                                                         %%%
%%%     E-mail: m.saadeghii@gmail.com                                                                      %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if nargin < 3
    error('Not enough input parameters!');
end

% parse input parameters

if (isfield(opts,'c'))
    c=opts.c;
else
    c=0.9;
end

if (isfield(opts,'alpha0'))
    alpha0=opts.alpha0;
else
    alpha0=500;
end

if (isfield(opts,'wF'))
    wF=opts.wF;
else
    wF=0.85;
end

if (isfield(opts,'wQ'))
    wQ=opts.wQ;
else
    wQ=0.85;
end

if (isfield(opts,'Ni'))
    Ni=opts.Ni;
else
    Ni=15;
end

if (isfield(opts,'T'))
    T=opts.T;
else
    T=3;
end

F=F0;   

N=size(F,2);
I=eye(N);
mu=zeros(1,No);

Z0=F'*F-I;
ZZ=abs(Z0);
alph=alpha0*max(abs(ZZ(:)));
Fo=F;

for iter=1:No
          
    mu(iter)=mucoh(F);

    for l=1:Ni
        
        % Q-update
        
        Z0o=Z0;
        Z2=(F'*F-I);
        Z0=Z2+wQ*(Z2-Z0o);
        
        z0=reshape(Z0,[N^2,1]);
        z1=ProjectOntoL1Ball(z0,alph);
        Z1=reshape(z1,[N,N]);
        
        Q=Z0-Z1+I;
        
        % F-update
        
        for j=1:T
            E=normc(F-muF*F*(F'*F-Q))+wF*(F-Fo);
            Fo=F;
            F=E;
        end
        
    end
        
    alph=alph*c;
    
end
