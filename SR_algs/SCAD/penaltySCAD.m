function [hbeta,hbetaprime]=penaltySCAD(beta,lambda,a)

%a=2;

if a <2
    error(' a should be >= 2')
end;
%lambda=1;

%
hbeta= lambda.*abs(beta).*(abs(beta)<=lambda);
hbeta= hbeta - (abs(beta).^2-2*a*lambda.*abs(beta)+ lambda.^2)/2/(a-1).*( (abs(beta)> lambda) & (abs(beta) <= a*lambda));
hbeta = hbeta + (a+1)*lambda.^2/2.* (abs(beta)> a*lambda);
%hbeta= abs(beta)-hbeta./lambda;


hbetaprime= lambda.*sign(beta) -((sign(beta).*lambda.*(abs(beta)<=lambda) + sign(beta).*lambda.*max(0,a*lambda-abs(beta))/(a-1)./lambda.*(abs(beta)>lambda)));