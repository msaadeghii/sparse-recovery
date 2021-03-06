function x=mycosamp(y,D,G,Dty,s,tr,maxiter)
m=size(D,2);
x=zeros(m,1);
x1=zeros(m,1);
iter=0;
ep=y'*y;
dep=0;
while ep > tr && iter<maxiter
   Gx=G*x;
   c=Dty-Gx; 
   [~, Gam]=selectsl(c,2*s);
   T=union(find(x),Gam);
   DT=D(:,T);
   DtD=DT'*DT;
   DTy=DT'*y;
   x1(T)= conjgrad(DtD,DTy,x(T));
   x1(~T)=0;
   [x ,~]=selectsl(x1,s);
   de=x'*Gx;
   ep=ep-de+dep;
   dep=de;
   iter=iter+1;
end

    function [yo ind]=selectsl(y,s)
        yab=abs(y);
        [a b]=sort(yab,'descend');
        yo=y;
        yo(yab<a(s))=0;
        ind=b(1:s);
       
    end
%     function z=mypinv()
%         M=DtD-eye(size(DtD));
%         z=x(T);
%        for i=1:3
%            z=DT'*y-M*z;
%        end
%     end
end