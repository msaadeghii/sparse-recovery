function mse=evalObjectiveVal(X,beta,y,lambda,penalty,param1)

datafit=0.5*(y-X*beta)'*(y-X*beta);
switch penalty
    case 'Log'
              mse=datafit+lambda*sum(log(abs(beta) +param1)) - lambda*log(param1);
    case 'Scad'
        scad=penaltySCAD(beta,lambda,param1);
        mse=datafit + sum( scad );
    case 'Zhang'
        
        mse=datafit + lambda*sum(penaltyZhang(beta,param1));
    case 'Lq'
        q=param1(1);
        epsilon=param1(2);
        mse =datafit + lambda * sum((abs(beta)+epsilon).^q);
end;