clear all
close all
clc

% Test CoSaMP sparse approximation algorithm
m =40;
n =20;
s =3;

scc1=0;
scc2=0;
scc3=0;

iter=100;
ti1=0;
ti2=0;
ti3=0;

trr=0.0001;
maxiter=100;

for i=1:iter
D = randn(n,m);
D=D*diag(sqrt(1./sum(D.*D)));
xo = zeros(m,1);
xo(randsample(m,s))= randn(s,1);
sn=0.01;
y=D*xo+random('norm',0,sn,[n 1]);
nxo=norm(xo);
D_pinv=pinv(D);

Dty=D'*y;
G=D'*D;
tr=1e-4;

mu_0=2;
sigma_min=0.00001;
sigma_decrease_factor=0.6;
L=5;

% CoSaMP
tic;x1=mycosamp(y,D,G,Dty,s,tr,maxiter);t1=toc;
% SP
tic;x2=mySP(y,D,G,Dty,s,tr,maxiter);t2=toc;
% SL0
tic;x3=SL0(D, y, sigma_min, sigma_decrease_factor, mu_0, L, D_pinv);t3=toc;

sc1=norm(x1-xo,'inf')<trr;
sc2=norm(x2-xo,'inf')<trr;
sc3=norm(x3-xo,'inf')<trr;

scc1=scc1+sc1;
scc2=scc2+sc2;
scc3=scc3+sc3;

ti1=ti1+t1;
ti2=ti2+t2;
ti3=ti3+t3;

display(i);
end

scr1=scc1/iter;
scr2=scc2/iter;
scr3=scc3/iter;

tim1=ti1/iter;
tim2=ti2/iter;
tim3=ti3/iter;

subplot(311);stem(xo,'r','filled');hold on;stem(x1,'b');legend('Original','Estimated');title(['CoSaMP ---- SRR: ',num2str(sc1),'--- Time: ',num2str(tim1)]);
subplot(312);stem(xo,'r','filled');hold on;stem(x2,'b');title(['SP ---- SRR: ',num2str(sc2),'--- Time: ',num2str(tim2)]);
subplot(313);stem(xo,'r','filled');hold on;stem(x3,'b');title(['SL0 ---- SRR: ',num2str(sc2),'--- Time: ',num2str(tim3)]);