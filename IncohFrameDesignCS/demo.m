% A demo of how to use the "IFD_AMPM" function for designing incoherent frames
close all;clear;clc

%% Experiment setup

n=15;   % dimension of frame vectors
m=50;  % number of frame vectors

mu_WB=sqrt((m-n)/(n*(m-1))); % Welch bound

F0=UNTF_iterative_projection(n, m, 50); % initialization with a UNTF

% IFD-AMPM parameters

muF=0.05;                  % step-size of the F-update problem
No=100;                      %  number of iterations
opts.c=0.9;               % threshold decaying factor (default: c=0.9)
opts.alpha0=500;      % factor of the initial value for alpha (default: alpha0=500)
opts.wF=0.85;            % weighting constant of F-update (default: wF=0.85)
opts.wQ=0.85;           % weighting constant of Q-update (default: wQ=0.85)
opts.Ni=15;                 % number of inner-loop iterations (default: Ni=15)
opts.T=3;                     % number of F-update iterations (default: T=3)

%% Algorithm

% IFD-AMPM
tic;[F, mu]=IFD_AMPM(F0, muF, No, opts);t_ifd=toc;

%% Plot the result

figure;plot(mu,'b','linewidth',2.5);hold on;
plot(1:No,mu_WB*ones(1,No),'--r','linewidth',2);
leg=legend('IFD-AMPM','WB');set(leg,'fontsize',15);
set(gca,'fontsize',15);
xlabel('iteration','FontName','Times New Roman','fontsize',15);
ylabel('$\mu(F)$','interpreter','latex','fontsize',15);
grid on;