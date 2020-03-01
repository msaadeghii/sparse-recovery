function [x,er,iter]=IPP_MCP(A, y, c, maxiter, Tf, L, A_pinv, xo,sigC,w,eps,gamma,noise_mode,rob_mode,a)

% This function implements the ISP-Hard algorithm proposed in the following
% paper:
% F. Ghayem, M. Sadeghi, M. Babaie-Zadeh, Ch. Jutten. "SL0 versus Iterative Hard/Soft Thresholding 
% for Sparse Signal Recovery: A Detailed Comparison and Accelerated
% Extensions", IEEE transaction on signal processing, 2016.
%
% for recovery of sparse signals from their compressed linear measurements:
% y=Ax+e: 
% x: the unknown sparse signal
% A: the measurement matrix
% y: the available measurements
% e: the noise vector
%
% Input parameters:
%                                   params.A :   Measurment matrix
%                                   params.y :    Measurnments
%                                   params.c:    Threshold decaying factor
%                                   params.maxiter:     Number of iterations as the iteration stop condition for threshold decaying loop 
%                                   params.Tf :    Minimum threshold level as the iteration stop condition for threshold decaying loop 
%                                   params.L :      Iteration number for sparsification_projection loop
%                                   params.A_pinv :      Pseudo_Inverse of measurment matrix
%                                   params.xo :     The unknown sparse signal
%                                   params.sigC :   Thresholding Initialization factor
%                                   params.w :  Extrapolation weight
%                                   params.eps :     Noise power
%                                   params.gamma :      Projection parameter using ADMM in robust case
%                                   params.noise_mode :     If equal to 1, then noisy case is considered
%                                   params.rob_mode :     If equal to 1, then robust case is considered
%                                   params.a :    SCAD parameter

%
% Output parameters:
%                                   params.x :      Recovered sparse signal
%                                   params.er :     Sparse signal recunstruction error
%                                   params.iter:    Number of passed thresholding iterations



x = A_pinv*y;            % Sparse signal initialization solving norm2 problem
xk=x;
xko=x;

% Main Loop

thr=sigC*max(abs(x));
er=[];
iter=1;

K=size(A,2);
M=inv(eye(K)+gamma*(A'*A));
lam=zeros(size(y));
epsp2=eps^2;

if noise_mode==0
    
    while iter <= maxiter || thr>=Tf  %thr > tol
        
        thr=c*thr;
        
        for l=1:L
            
            yk=xk+w*(xk-xko);                            % Extrapolation
            x=MCP_Prox(yk,thr,a);       % Thresholding
            x = x - A_pinv*(A*x-y);                      % Projection
            xko=xk;
            xk=x;
            
        end
        
        rel_res=norm(x-xo)/norm(xo);
        er(iter)=rel_res;
        iter=iter+1;
        
    end
    
elseif rob_mode==0
    
    while iter <= maxiter || thr>=Tf  %thr > tol
        
        thr=c*thr;
        
        for l=1:L
            
            yk=xk+w*(xk-xko);                            % Extrapolation
            x=MCP_Prox(yk,thr,a);                 % Thresholding
            zx=y-A*x;
            
            if (zx'*zx) > epsp2
                x = x - A_pinv*(A*x-y);                  % Eftekhari's Robust Projection
            end
            
            xko=xk;
            xk=x;
            
        end
        
        rel_res=norm(x-xo)/norm(xo);
        er(iter)=rel_res;
        iter=iter+1;
        
    end
    
elseif rob_mode==1
    
    while iter <= maxiter || thr>=Tf  %thr > tol
        
        thr=c*thr;
        
        for l=1:L
            
            yk=xk+w*(xk-xko);                            % Extrapolation
            x=MCP_Prox(yk,thr,a);       % Thresholding
            zx=y-A*x;
            x = x0;
			
            while (zx'*zx) > epsp2                 % Proposed Robust Projection
                
                % z-update
                z=y-A*x+1/gamma*lam;
                zn=sqrt(z'*z);
                
                if zn > eps
                    z=z/zn *eps;                         % z projection to A_z set
                end
                
                % x-update
                G=y-z+1/gamma*lam;
                x=M*(x0+gamma*(A'*G));
                
                zdif=z-y+A*x;
                
                % lam-update
                lam=lam-gamma*(zdif);
                zx=z-zdif;
                
            end
            
            xko=xk;
            xk=x;
            
        end
        
        rel_res=norm(x-xo)/norm(xo);
        er(iter)=rel_res;
        iter=iter+1;
        
    end
    
    
end

end