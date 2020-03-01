function [BetaCCCP,alphaAdapA,BetaHist,t,Mse]=DCLasso(Kapp,yapp,lambda0,param1,type,CCCPoption,Lassooption)

 

% CCCP lasso - reweigthed lasso using
%
%
tic;

TailleDic=size(Kapp,2);
% lambda = lambda0;
if ~isfield(CCCPoption,'init');
    BetaCCCP = zeros(TailleDic,1);
    lambda=lambda0;%*ones(TailleDic,1);
else
    BetaCCCP = CCCPoption.init;
    lambda = weights(abs(BetaCCCP),lambda0,param1,type);
end;

BetaOld=BetaCCCP;BetaOld(1)=1;
BetaHist=[];
nbiter=0;
alphaAdapA=0;

while norm(BetaCCCP-BetaOld,'inf')>CCCPoption.CCCPstopcrit & nbiter < CCCPoption.nbitermax
    nbiter=nbiter+1;
    BetaOld=BetaCCCP;
    
    switch CCCPoption.algo
        case 'coord'
    [BetaCCCP,alphaAdapA]=sparseleastsquare(Kapp,yapp,lambda,Lassooption);
        case 'GP'
    [BetaCCCP,alphaAdapA]=LassoGP(Kapp,yapp,lambda,Lassooption);
    end;
    lambda=weights(abs(BetaCCCP),lambda0,param1,type);
    TimeHist(nbiter)=toc;
    Lassooption.alphainit=BetaCCCP;
    BetaHist=[BetaHist BetaCCCP];
    Mse(nbiter)=evalObjectiveVal(Kapp,BetaCCCP,yapp,lambda0,type,param1);
end;

t = toc;


%--------------------------------------------------------------------------
%
%               Others functions
%
%--------------------------------------------------------------------------


%---------------------------------------------------------
%           Process the weights at each iteration
%---------------------------------------------------------
function  lambda=weights(vec,lambda0,param1,type);
switch type
    case 'Zhang'
        [hbeta,hbetaprime]=penaltyZhang(vec,param1);
        lambda=[lambda0.*(1-hbetaprime)];%;lambda0.*(1+hbetaprime)];
    case 'Scad'
        [hbeta,hbetaprime]=penaltySCAD(vec,lambda0,param1);
        hbetaprime=hbetaprime./lambda0;
        lambda=[lambda0.*(1-hbetaprime)];%;lambda0.*(1+hbetaprime)];
    case 'Log'
       lambda=[lambda0./(vec+param1)];%;lambda0./(abs(BetaCCCP)+param1)];      
    case 'Lq'
        if length(param1)==1
        q=param1;
        epsilon=0;
        else 
            q=param1(1);
            epsilon=param1(2);
        end;
        lambda=[lambda0.*q*(vec + epsilon).^(q-1)];%;lambda0.*q*(abs(BetaCCCP) + epsilon).^(q-1)];
end,


