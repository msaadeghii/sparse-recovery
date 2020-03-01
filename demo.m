%% This script compares the performance of different sparse recovery algorithms from noisy measurements.

clear all;close all;clc

addpath(genpath(cd))

tic;

% Sensing matrix dimensions n*m
m = 1000;
n = 400;

k = [50 75 100:20:140 150:10:200 205:5:250]; %Sparsity Level

SigmaOn = 1; % Non-zeros standard deviation
SigmaNoise = 0.01; % Noise standard deviation
path(path, 'Other Solvers');

% Number of Monte-Carlo simulations
NMont = 500;

%SCSA Params
eps1_SCSA = 1e-3;
eps1_min_SCSA = 1e-4;
eps2_SCSA = 1e-2;
eps2_min_SCSA = 1e-3;
c_SCSA = 0.1;

% ISP_SCAD Params
maxiter=300;             % Number of iterations as the iteration stop condition for threshold decaying loop
Tf=1e-15;                   % Minimum threshold level as the iteration stop condition for threshold decaying loop
L=3;                            % Iteration number for sparsification_projection loop
c=0.9;                         % Threshold decaying factor
sigC=3;                      % Thresholding Initialization factor
a=30;                          % SCAD parameter
ww=0.95;                                              % Extrapolation weight
gamma=0.4;                                        % Projection parameter using ADMM in robust case
noise_mode=1;                                % Noise mode parameter
rob_mode=1;                                    % robust mode parameter

lambda_BPDN = SigmaNoise * 2 * 1.05 * norminv(1 - 0.5 / (2 * m),0,1);

% SCAD Params
a_SCAD = 3.7;
lambda_SCAD = 0.5 * lambda_BPDN;
eps_min_SCAD = 1e-4;
eps_SCAD = 1e-3;
CCCPoption.CCCPstopcrit = 1e-4;
CCCPoption.nbitermax = 20;
CCCPoption.algo = 'GP';
Lassooption.bias = 0;
Lassooption.verbose = 0;
Lassooption.KKTcond = 0.001;
Lassooption.nbitermax = 10000;

% EM-GM-AMP Params
optEM.heavy_tailed = false;  % since operating on a sparse signal
optGAMP.nit = 20;
optGAMP.tol = 1e-5;
optEM.maxTol = 1e-5;
optEM.maxEMiter = 20;

% Initialize matrices
SNR_IPP_SCAD = ones(length(k),NMont);
SNR_EMGMAMP = ones(length(k),NMont);
SNR_SCSA_FIT = ones(length(k),NMont);
SNR_GOMP = ones(length(k),NMont);
SNR_SCAD = ones(length(k),NMont);

SE_IPP_SCAD = ones(length(k),NMont);
SE_EMGMAMP = ones(length(k),NMont);
SE_SCSA_FIT = ones(length(k),NMont);
SE_GOMP = ones(length(k),NMont);
SE_SCAD = ones(length(k),NMont);

SupRate_ISP_SCAD = zeros(length(k),1);
SupRate_EMGMAMP = zeros(length(k),1);
SupRate_SCSA_FIT = zeros(length(k),1);
SupRate_GOMP = zeros(length(k),1);
SupRate_SCAD = zeros(length(k),1);

Car_IPP_SCAD = zeros(length(k),NMont);
Car_EMGMAMP = zeros(length(k),NMont);
Car_SCSA_FIT = zeros(length(k),NMont);
Car_GOMP = zeros(length(k),NMont);
Car_SCAD = zeros(length(k),NMont);

Time_IPP_SCAD = zeros(length(k),NMont);
Time_EMGMAMP= zeros(length(k),NMont);
Time_SCSA_FIT = zeros(length(k),NMont);
Time_GOMP = zeros(length(k),NMont);
Time_SCAD = zeros(length(k),NMont);

SNR_Oracle = inf * ones(length(k),NMont);
SE_Oracle = inf * ones(length(k),NMont);

parpool(2);

for i = 1 : length(k)
    
    txt = sprintf('********************* Simulation for k = %d ************************',k(i));
    disp(txt);
    
    parfor j = 1 : NMont
        
        % Generate sparse signal
        xber = zeros(m,1);
        supp = randperm(m, k(i))';
        xber(supp) = 1;        
        x1 = normrnd(0, SigmaOn, m, 1);
        x  = xber .* x1;
        
        % Generate measurement matrix
        A = normrnd(0, 1, n, m);
        W=zeros(m);
        for p=1:m
            W(p,p)=1/norm(A(:,p));
        end
        
        Ah=A*W;
        
        % Generate measurements
        w = normrnd(0,SigmaNoise,n,1);
        y = A * x + w;
        
        A_pinv=pinv(Ah);         % Pseudo_Inverse of random measurment matrix
        tol=norm(w);                         	     % Random noise power    
        
        %% IPP-SCAD
        tic;  [xhat] =IPP_SCAD(Ah, y, c, maxiter, Tf, L, A_pinv, x, sigC, ww, tol, gamma, noise_mode, rob_mode, a); t_IPP_SCAD=toc; Time_IPP_SCAD(i,j)=t_IPP_SCAD;
        xhat=W*xhat;
        SNR_IPP_SCAD(i,j) = -20*log10( norm(xhat - x,2) / norm(x,2) );
        SE_IPP_SCAD(i,j) = norm(xhat - x,2)^2;
        [~,indx1] = sort(abs(xhat),'descend');
        indx1 = indx1(1:k(i));
        Car_IPP_SCAD(i,j) = length( find( sort(indx1) == sort(supp) ) ) / k(i);
        
        %% EM-GM-AMP
        tic;[xhat, EMfin] = EMGMAMP(y, Ah, optEM,optGAMP);t_EMGM=toc;Time_EMGMAMP(i,j)=t_EMGM;
        xhat=W*xhat;
        SNR_EMGMAMP(i,j) = -20*log10( norm(xhat - x,2) / norm(x,2) );
        SE_EMGMAMP(i,j) = norm(xhat - x,2)^2;
        [~,indx1] = sort(abs(xhat),'descend');
        indx1 = indx1(1:k(i));
        Car_EMGMAMP(i,j) = length( find( sort(indx1) == sort(supp) ) ) / k(i);
        
        %% SCSA-FIT
        [xhat,Time_SCSA_FIT(i,j)] = SCSA_FIT(Ah, y, min(eps1_min_SCSA, eps1_SCSA * lambda_BPDN), 0.1 * min(eps2_min_SCSA, eps2_SCSA * lambda_BPDN), c_SCSA,lambda_BPDN);
        xhat=W*xhat;
        SNR_SCSA_FIT(i,j) = -20*log10( norm(xhat - x,2) / norm(x,2) );
        SE_SCSA_FIT(i,j) = norm(xhat - x,2)^2;
        [~,indx1] = sort(abs(xhat),'descend');
        indx1 = indx1(1:k(i));
        Car_SCSA_FIT(i,j) = length( find( sort(indx1) == sort(supp) ) ) / k(i);
        
        %% GOMP     
        tic;[xhat] = islsp_EstgOMP(y, Ah, 1000, 4, tol);t_GOMP=toc;Time_GOMP(i,j)=t_GOMP;
        xhat=W*xhat;
        SNR_GOMP(i,j) = -20*log10( norm(xhat - x,2) / norm(x,2) );
        SE_GOMP(i,j) = norm(xhat - x,2)^2;
        [~,indx1] = sort(abs(xhat),'descend');
        indx1 = indx1(1:k(i));
        Car_GOMP(i,j) = length( find( sort(indx1) == sort(supp) ) ) / k(i);
        
        
        %% SCAD     
         [xhat,~,~,Time_SCAD(i,j),~] = DCLasso(Ah, y, lambda_SCAD, a_SCAD, 'Scad', CCCPoption, Lassooption);
        xhat=W*xhat;
         SNR_SCAD(i,j) = -20*log10( norm(xhat - x,2) / norm(x,2) );
        SE_SCAD(i,j) = norm(xhat - x,2)^2;
        indx2 = find( abs(xhat) > SigmaNoise);
        [~,indx1] = sort(abs(xhat),'descend');
        indx1 = indx1(1:k(i));
        Car_SCAD(i,j) = length( find( sort(indx1) == sort(supp) ) ) / k(i);
        
        if (rem(j,ceil(NMont/10)) == 0)
            txt = sprintf('Monte Carlo No. %d / %d.',j,NMont);
            disp(txt);
        end;
        
    end;
    
    % Save the results
%     save(['MSE_ALL_Algorithms_OurNormalization,','_','noise_level_',num2str(SigmaNoise),'.mat']...
%         ,'Time_ISP_SCAD','SNR_ISP_SCAD','SE_ISP_SCAD','Car_ISP_SCAD'...
%         ,'Time_SCSA_FIT','SNR_SCSA_FIT','SE_SCSA_FIT','Car_SCSA_FIT'...
%         ,'Time_SCAD','SNR_SCAD','SE_SCAD','Car_SCAD'...
%         ,'Time_EMGMAMP','SNR_EMGMAMP','Car_EMGMAMP'...
%         ,'Time_GOMP','SNR_GOMP','SE_GOMP','Car_GOMP');
    
end;

MSE_ISP_SCAD = mean(SE_IPP_SCAD,2)';
MSE_EMGMAMP = mean(SE_EMGMAMP,2)';
MSE_SCSA_FIT = mean(SE_SCSA_FIT,2)';
MSE_GOMP = mean(SE_GOMP,2)';
MSE_SCAD = mean(SE_SCAD,2)';

% delete(gcp('nocreate'))

%%
close all;

for i = 1 : length(k),
    SupRate_ISP_SCAD(i) = sum(Car_IPP_SCAD(i,:)) / NMont;
    SupRate_SCAD(i) = sum(Car_SCAD(i,:)) / NMont;
    SupRate_SCSA_FIT(i) = sum(Car_SCSA_FIT(i,:)) / NMont;
    SupRate_GOMP(i) = sum(Car_GOMP(i,:)) / NMont;
    SupRate_EMGMAMP(i) = sum(Car_EMGMAMP(i,:)) / NMont;
end;

figure,
plot(k,mean(SNR_IPP_SCAD,2),'r',k,mean(SNR_EMGMAMP,2),'.-m',k,mean(SNR_SCSA_FIT,2),'c',k,mean(SNR_GOMP,2),'b',k,mean(SNR_SCAD,2),'.-b');
legend('ISP-SCAD','EMGMAMP','SCSA-FIT','GOMP','SCAD');
title('MSNR');

figure,
plot(k,SupRate_ISP_SCAD,'--b',k,SupRate_EMGMAMP,'b',k,SupRate_SCSA_FIT,'--r',k,SupRate_GOMP,'r',...
    k,SupRate_SCAD,'c');
legend('ISP-SCAD','EMGMAMP','SCSA-FIT','GOMP','SCAD');
title('Support Recovery Rate');

figure,
semilogy(k,mean(Time_IPP_SCAD,2),'r',k,mean(Time_EMGMAMP,2),'.-r',k,mean(Time_SCSA_FIT,2),'b',k,mean(Time_GOMP,2),'.-b',...
    k,mean(Time_SCAD,2),'.-g');
legend('ISP-SCAD','EMGMAMP','SCSA-FIT','GOMP','SCAD');
title('Average Run Time');

t = toc
